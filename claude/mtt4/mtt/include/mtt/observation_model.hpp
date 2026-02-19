#pragma once
#include <algorithm>
#include "types.hpp"

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// Convert a RawMeasurement (degrees) into a filter-ready Measurement (radians)
// using the globally uniform noise covariance from MHTConfig.
// ─────────────────────────────────────────────────────────────────────────────
inline Measurement convert_measurement(const RawMeasurement& raw,
                                       double noise_deg2) noexcept {
    Measurement m;
    m.timestamp  = raw.timestamp;
    m.sensor_id  = raw.sensor_id;

    // Angles to radians
    m.z(0) = deg2rad(raw.az_deg);
    m.z(1) = deg2rad(raw.el_deg);

    // Covariance scales by square of conversion factor
    const double s = DEG2RAD * DEG2RAD;
    m.R = Mat2d::Identity() * (noise_deg2 * s);   // diagonal σ² diag
    return m;
}

// ─────────────────────────────────────────────────────────────────────────────
// Non-linear observation function h(state, sensor)
// Returns [azimuth_rad, elevation_rad] predicted from the 6D ECEF state.
// ─────────────────────────────────────────────────────────────────────────────
inline Vec2d bearing_observation(const Vec6d& state,
                                 const SensorInfo& sensor) noexcept {
    // Line-of-sight vector from sensor to target in ECEF
    const Vec3d r = state.head<3>() - sensor.ecef;
    const double range = r.norm();
    if (range < 1.0) {
        // Target at sensor location – undefined; return zeros gracefully
        return Vec2d::Zero();
    }

    // Project onto ENU axes
    const double e_proj = sensor.east_ecef.dot(r);
    const double n_proj = sensor.north_ecef.dot(r);
    const double u_proj = sensor.up_ecef.dot(r);

    const double az  = std::atan2(e_proj, n_proj);  // [-π, π]
    const double el  = std::asin( std::clamp(u_proj / range, -1.0, 1.0));
    return { az, el };
}

// ─────────────────────────────────────────────────────────────────────────────
// Innovation (wrapped azimuth difference to handle ±π boundary)
// ─────────────────────────────────────────────────────────────────────────────
inline Vec2d compute_innovation(const Vec2d& z_measured,
                                const Vec2d& z_predicted) noexcept {
    Vec2d nu = z_measured - z_predicted;
    // Wrap azimuth component into [-π, π]
    while (nu(0) >  M_PI) nu(0) -= 2.0*M_PI;
    while (nu(0) < -M_PI) nu(0) += 2.0*M_PI;
    return nu;
}

// ─────────────────────────────────────────────────────────────────────────────
// Analytical Jacobian H = ∂h/∂x  (2×6)
// Included for reference / EKF fall-back; not used in the primary UKF path.
// Caution: undefined near zenith (el → ±π/2).
// ─────────────────────────────────────────────────────────────────────────────
inline Mat2x6 bearing_jacobian(const Vec6d& state,
                               const SensorInfo& sensor) noexcept {
    const Vec3d r    = state.head<3>() - sensor.ecef;
    const double e_p = sensor.east_ecef.dot(r);
    const double n_p = sensor.north_ecef.dot(r);
    const double u_p = sensor.up_ecef.dot(r);
    const double rho2= r.squaredNorm();
    const double rho_h2 = e_p*e_p + n_p*n_p;   // horizontal range²
    const double rho_h  = std::sqrt(rho_h2);

    Mat2x6 H = Mat2x6::Zero();

    // ∂az/∂pos
    for (int i = 0; i < 3; ++i) {
        const double de_i = sensor.east_ecef(i);
        const double dn_i = sensor.north_ecef(i);
        H(0, i) = (rho_h2 > 1e-6) ?
                  (de_i * n_p - dn_i * e_p) / rho_h2 : 0.0;
    }

    // ∂el/∂pos
    for (int i = 0; i < 3; ++i) {
        const double du_i = sensor.up_ecef(i);
        const double dterm = (rho2 > 1e-6) ?
                             (du_i * rho2 - u_p * r(i)) / (rho2 * rho_h) : 0.0;
        H(1, i) = (rho_h > 1e-6) ? dterm : 0.0;
    }
    // Velocity columns remain zero (bearing-only)
    return H;
}

}  // namespace mtt
