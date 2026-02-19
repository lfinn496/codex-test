#pragma once
#include "types.hpp"
#include "observation_model.hpp"
#include <Eigen/Cholesky>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// UKF state structure (owned by each Track object)
// ─────────────────────────────────────────────────────────────────────────────
struct UKFState {
    Vec6d  x = Vec6d::Zero();         // mean state [px,py,pz,vx,vy,vz]
    Mat6d  P = Mat6d::Identity();     // covariance
    double timestamp = 0.0;           // time of last update [s]
};

// ─────────────────────────────────────────────────────────────────────────────
// Unscented Kalman Filter (Merwe scaled sigma-point variant)
// ─────────────────────────────────────────────────────────────────────────────
class UKF {
public:
    static constexpr int N  = 6;                // state dimension
    static constexpr int NS = 2*N + 1;          // number of sigma points (13)

    // Sigma-point scaling parameters
    double alpha  = 1e-3;
    double beta   = 2.0;
    double kappa  = 0.0;

    // CWNA process noise PSD [m²/s³]
    double psd    = 3.0;

    // Constructor applies config
    explicit UKF(const MHTConfig& cfg)
        : alpha(cfg.ukf_alpha), beta(cfg.ukf_beta), kappa(cfg.ukf_kappa),
          psd(cfg.psd_default) {
        compute_weights();
    }

    // ─── CWNA Prediction ─────────────────────────────────────────────────────
    // Propagates x and P forward by dt seconds.
    void predict(UKFState& s, double dt) const {
        // State transition: constant-velocity (CWNA)
        const Mat6d F = build_F(dt);
        const Mat6d Q = build_Q(dt, psd);
        // Standard linear propagation for CWNA (non-linear not required here)
        s.x = F * s.x;
        s.P = F * s.P * F.transpose() + Q;
        // Ensure P stays symmetric positive-definite
        s.P = 0.5*(s.P + s.P.transpose());
    }

    // ─── Sequential UKF Update ───────────────────────────────────────────────
    // Updates state given one bearing measurement from one sensor.
    // Returns the normalised innovation (Mahalanobis) squared for gating/scoring.
    double update(UKFState& s, const Measurement& meas,
                  const SensorInfo& sensor) const {
        // 1. Generate sigma points
        Eigen::Matrix<double, N, NS> sigmas;
        generate_sigma_points(s, sigmas);

        // 2. Propagate sigma points through h()
        Eigen::Matrix<double, 2, NS> z_sigmas;
        for (int i = 0; i < NS; ++i) {
            Vec6d xi;
            xi = sigmas.col(i);
            z_sigmas.col(i) = bearing_observation(xi, sensor);
        }

        // 3. Predicted measurement mean
        Vec2d z_hat = Vec2d::Zero();
        for (int i = 0; i < NS; ++i) z_hat += Wm_(i) * z_sigmas.col(i);
        // Wrap azimuth
        while (z_hat(0) >  M_PI) z_hat(0) -= 2.0*M_PI;
        while (z_hat(0) < -M_PI) z_hat(0) += 2.0*M_PI;

        // 4. Innovation covariance S and cross-covariance Pxz
        Mat2d  S  = meas.R;
        Mat6x2 Pxz = Mat6x2::Zero();
        for (int i = 0; i < NS; ++i) {
            Vec2d dz = compute_innovation(z_sigmas.col(i), z_hat);
            Vec6d dx = sigmas.col(i) - s.x;
            S   += Wc_(i) * (dz * dz.transpose());
            Pxz += Wc_(i) * (dx * dz.transpose());
        }
        S = 0.5*(S + S.transpose());

        // 5. Innovation
        const Vec2d nu = compute_innovation(meas.z, z_hat);

        // 6. Mahalanobis distance² for gating
        const Mat2d S_inv = S.inverse();
        const double maha2 = nu.transpose() * S_inv * nu;

        // 7. Kalman gain and state update
        const Mat6x2 K = Pxz * S_inv;
        s.x = s.x + K * nu;
        s.P = s.P - K * S * K.transpose();
        s.P = 0.5*(s.P + s.P.transpose());

        return maha2;
    }

    // Compute the innovation covariance and Mahalanobis distance WITHOUT
    // modifying state (for gating).
    double gating_distance(const UKFState& s, const Measurement& meas,
                           const SensorInfo& sensor) const {
        Eigen::Matrix<double, N, NS> sigmas;
        generate_sigma_points(s, sigmas);

        Eigen::Matrix<double, 2, NS> z_sigmas;
        for (int i = 0; i < NS; ++i) {
            Vec6d xi; xi = sigmas.col(i);
            z_sigmas.col(i) = bearing_observation(xi, sensor);
        }
        Vec2d z_hat = Vec2d::Zero();
        for (int i = 0; i < NS; ++i) z_hat += Wm_(i) * z_sigmas.col(i);

        Mat2d S = meas.R;
        for (int i = 0; i < NS; ++i) {
            Vec2d dz = compute_innovation(z_sigmas.col(i), z_hat);
            S += Wc_(i) * (dz * dz.transpose());
        }
        S = 0.5*(S + S.transpose());
        const Vec2d nu = compute_innovation(meas.z, z_hat);
        return nu.transpose() * S.inverse() * nu;
    }

    // ─── Build CWNA matrices ──────────────────────────────────────────────────
    static Mat6d build_F(double dt) noexcept {
        Mat6d F = Mat6d::Identity();
        F.topRightCorner<3,3>() = Mat3d::Identity() * dt;
        return F;
    }

    static Mat6d build_Q(double dt, double q) noexcept {
        // Continuous White Noise Acceleration discretised
        const double dt2 = dt*dt;
        const double dt3 = dt2*dt;
        Mat6d Q = Mat6d::Zero();
        Q.topLeftCorner<3,3>()     = Mat3d::Identity() * (q * dt3 / 3.0);
        Q.topRightCorner<3,3>()    = Mat3d::Identity() * (q * dt2 / 2.0);
        Q.bottomLeftCorner<3,3>()  = Mat3d::Identity() * (q * dt2 / 2.0);
        Q.bottomRightCorner<3,3>() = Mat3d::Identity() * (q * dt);
        return Q;
    }

private:
    // Pre-computed weights
    Eigen::Matrix<double, NS, 1> Wm_, Wc_;

    void compute_weights() {
        const double lambda = alpha*alpha*(N + kappa) - N;
        Wm_(0) = lambda / (N + lambda);
        Wc_(0) = Wm_(0) + (1.0 - alpha*alpha + beta);
        const double w = 1.0 / (2.0*(N + lambda));
        for (int i = 1; i < NS; ++i) { Wm_(i) = w; Wc_(i) = w; }
    }

    void generate_sigma_points(const UKFState& s,
                               Eigen::Matrix<double, N, NS>& sigmas) const {
        const double lambda = alpha*alpha*(N + kappa) - N;
        // Cholesky of (N+λ)P
        Eigen::LLT<Mat6d> llt(( N + lambda) * s.P);
        const Mat6d L = llt.matrixL();

        sigmas.col(0) = s.x;
        for (int i = 0; i < N; ++i) {
            sigmas.col(i+1)   = s.x + L.col(i);
            sigmas.col(i+1+N) = s.x - L.col(i);
        }
    }
};

}  // namespace mtt
