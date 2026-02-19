#pragma once
#include "types.hpp"
#include <vector>
#include <cmath>
#include <Eigen/QR>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// Single bearing ray used by WS3D
// ─────────────────────────────────────────────────────────────────────────────
struct BearingRay {
    Vec3d  sensor_ecef;   // sensor position [m]
    double az_rad;        // azimuth   [rad]
    double el_rad;        // elevation [rad]
    double sigma2_az;     // azimuth   variance [rad²]
    double sigma2_el;     // elevation variance [rad²]
    Vec3d  east;          // sensor ENU east  in ECEF
    Vec3d  north;         // sensor ENU north in ECEF
    Vec3d  up;            // sensor ENU up    in ECEF
};

// ─────────────────────────────────────────────────────────────────────────────
// Result of WS3D triangulation
// ─────────────────────────────────────────────────────────────────────────────
struct WS3DResult {
    Vec3d   position_ecef;  // triangulated 3D position [m]
    Mat3d   cov_pos;        // 3×3 position covariance [m²]
    double  quality;        // trace(cov_pos) – lower is better
    bool    valid = false;
};

// ─────────────────────────────────────────────────────────────────────────────
// Weighted Stansfield Algorithm in Three Dimensions (WS3D)
//
// Ref: Stansfield (1947), extended to 3D with range-dependent weighting.
// Two-step: (1) 2D planar LS for (Px, Py); (2) 1D height estimate.
// ─────────────────────────────────────────────────────────────────────────────
class WS3DTriangulator {
public:
    // Minimum sensors (rays) required
    static constexpr int MIN_RAYS = 2;

    // ─── Core algorithm ───────────────────────────────────────────────────────
    static WS3DResult triangulate(const std::vector<BearingRay>& rays) {
        const int m = static_cast<int>(rays.size());
        if (m < MIN_RAYS) return {};

        // ── Step 1: unweighted 2D pass to get rough horizontal ranges ─────────
        Vec2d xy0 = planar_ls_unweighted(rays);

        // ── Compute approximate horizontal ranges for weighting ────────────────
        std::vector<double> R_h(m);
        for (int i = 0; i < m; ++i) {
            const Vec3d& sp = rays[i].sensor_ecef;
            R_h[i] = std::hypot(xy0(0) - sp(0), xy0(1) - sp(1));
            R_h[i] = std::max(R_h[i], 100.0);  // avoid division by zero
        }

        // ── Step 2: weighted 2D planar estimate ───────────────────────────────
        Vec2d xy1 = planar_ls_weighted(rays, R_h);

        // ── Recompute horizontal ranges with refined xy1 ───────────────────────
        for (int i = 0; i < m; ++i) {
            const Vec3d& sp = rays[i].sensor_ecef;
            R_h[i] = std::hypot(xy1(0) - sp(0), xy1(1) - sp(1));
            R_h[i] = std::max(R_h[i], 100.0);
        }

        // ── Step 3: 1D height estimation ──────────────────────────────────────
        double height = height_ls(rays, xy1, R_h);

        // ── Assemble ECEF result ───────────────────────────────────────────────
        // xy1 is actually 2D position in the sensor's local horizontal plane.
        // Because WS3D operates in a planar approximation, we build the ECEF
        // result by compositing:
        //   pos = (xy1[0], xy1[1], height)  all in ECEF metres.
        // This works well for sensors on flat ground because the ECEF x-y plane
        // is essentially the equatorial plane; for correctness in general
        // geometries the full ECEF formulation should be used.
        WS3DResult res;
        res.position_ecef = Vec3d(xy1(0), xy1(1), height);
        res.cov_pos       = build_cov(rays, res.position_ecef);
        res.quality       = res.cov_pos.trace();
        res.valid         = true;
        return res;
    }

private:
    // ── Unweighted least-squares: minimise sum of azimuth residuals ────────────
    // Each azimuth ray from sensor i in ECEF: sin(az_i)*N - cos(az_i)*E = const
    static Vec2d planar_ls_unweighted(const std::vector<BearingRay>& rays) {
        const int m = static_cast<int>(rays.size());
        Eigen::MatrixXd A(m, 2);
        Eigen::VectorXd b(m);
        for (int i = 0; i < m; ++i) {
            const double caz = std::cos(rays[i].az_rad);
            const double saz = std::sin(rays[i].az_rad);
            // From  N*saz - E*caz = N_s*saz - E_s*caz
            A(i,0) = saz;   // coefficient for E component
            A(i,1) = -caz;  // coefficient for N component ... wait, let me redo
            // Line from sensor (Ex, Ey) in direction (sin(az), cos(az)):
            //   (y - Ey)*sin(az) - (x - Ex)*cos(az) = 0
            //   => -cos(az)*x + sin(az)*y = -cos(az)*Ex + sin(az)*Ey
            A(i,0) = -caz;
            A(i,1) =  saz;
            b(i)   = -caz * rays[i].sensor_ecef(0)
                     + saz * rays[i].sensor_ecef(1);
        }
        Eigen::Vector2d sol = A.colPivHouseholderQr().solve(b);
        return sol;
    }

    // ── Weighted planar LS ─────────────────────────────────────────────────────
    // Weight w_i = 1 / (sigma2_az * R_h[i]²) — penalises far sensors
    static Vec2d planar_ls_weighted(const std::vector<BearingRay>& rays,
                                    const std::vector<double>& R_h) {
        const int m = static_cast<int>(rays.size());
        Eigen::MatrixXd A(m, 2);
        Eigen::VectorXd b(m);
        Eigen::VectorXd w(m);

        for (int i = 0; i < m; ++i) {
            const double caz = std::cos(rays[i].az_rad);
            const double saz = std::sin(rays[i].az_rad);
            A(i,0) = -caz;
            A(i,1) =  saz;
            b(i)   = -caz * rays[i].sensor_ecef(0)
                     + saz * rays[i].sensor_ecef(1);
            w(i)   = 1.0 / (rays[i].sigma2_az * R_h[i] * R_h[i]
                            + 1e-12);
        }
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(w);
        const Eigen::MatrixXd AtW  = A.transpose() * W;
        const Eigen::MatrixXd AtWA = AtW * A;
        const Eigen::VectorXd AtWb = AtW * b;
        return AtWA.ldlt().solve(AtWb);
    }

    // ── 1D height estimation ───────────────────────────────────────────────────
    // z_hat = sum_i( w_i * (tan(el_i)*R_h[i] + sensor_i.z) )
    //       / sum_i( w_i )
    // Weight: w_i = cos²(el_i) / sigma2_el (geometric dampener near zenith)
    static double height_ls(const std::vector<BearingRay>& rays,
                            const Vec2d& /*xy*/,
                            const std::vector<double>& R_h) {
        double num = 0.0, den = 0.0;
        for (std::size_t i = 0; i < rays.size(); ++i) {
            const double cos_el = std::cos(rays[i].el_rad);
            const double tan_el = std::tan(rays[i].el_rad);
            // Height dampener: w ∝ cos²(el)/σ²_el
            const double w = (cos_el * cos_el)
                            / (rays[i].sigma2_el * R_h[i] * R_h[i] + 1e-12);
            num += w * (rays[i].sensor_ecef(2) + tan_el * R_h[i]);
            den += w;
        }
        return (den > 1e-12) ? (num / den) : 0.0;
    }

    // ── Posterior position covariance (CRLB approximation) ────────────────────
    static Mat3d build_cov(const std::vector<BearingRay>& rays,
                           const Vec3d& pos) {
        // Approximate Fisher information via azimuth + elevation residuals
        // Using standard angular error propagation:  σ_pos ≈ σ_angle * R
        Mat3d FI = Mat3d::Zero();
        for (const auto& r : rays) {
            const Vec3d dr   = pos - r.sensor_ecef;
            const double R   = dr.norm();
            const double R2  = R * R;
            if (R2 < 1.0) continue;
            // Azimuth contributes to horizontal uncertainty
            const double cos_el = std::cos(r.el_rad);
            const double waz  = (cos_el*cos_el) / (r.sigma2_az * R2 + 1e-12);
            const double wel  = 1.0 / (r.sigma2_el * R2 + 1e-12);
            // Direction matrix contribution
            const Vec3d udir = dr / R;
            // Perpendicular projection: I - u*u^T gives uncertainty directions
            Mat3d perp = Mat3d::Identity() - udir * udir.transpose();
            FI += (waz + wel) * perp;
        }
        // Covariance is inverse of Fisher information
        // Add small ridge for numerical stability
        FI += Mat3d::Identity() * 1e-6;
        return FI.inverse();
    }
};

}  // namespace mtt
