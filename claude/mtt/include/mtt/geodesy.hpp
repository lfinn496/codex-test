#pragma once
#include "types.hpp"

namespace mtt {
namespace geodesy {

// ─────────────────────────────────────────────────────────────────────────────
// WGS-84: prime vertical radius of curvature N(φ)
// ─────────────────────────────────────────────────────────────────────────────
inline double prime_vertical_radius(double lat_rad) noexcept {
    const double sin_lat = std::sin(lat_rad);
    return wgs84::a / std::sqrt(1.0 - wgs84::e2 * sin_lat * sin_lat);
}

// ─────────────────────────────────────────────────────────────────────────────
// Geodetic (lat°, lon°, alt_m) → ECEF [m]
// All trig matches the standard ITRF/WGS-84 formulation.
// ─────────────────────────────────────────────────────────────────────────────
Vec3d lla_to_ecef(double lat_deg, double lon_deg, double alt_m) noexcept;

// ─────────────────────────────────────────────────────────────────────────────
// ECEF [m] → Geodetic (lat°, lon°, alt_m)
// Bowring/Heikkinen iterative solution – converges in 3 iterations.
// ─────────────────────────────────────────────────────────────────────────────
void ecef_to_lla(const Vec3d& ecef,
                 double& lat_deg, double& lon_deg, double& alt_m) noexcept;

// ─────────────────────────────────────────────────────────────────────────────
// Compute ENU unit vectors in ECEF frame for a sensor at (lat°, lon°)
// ─────────────────────────────────────────────────────────────────────────────
void compute_enu_axes(double lat_deg, double lon_deg,
                      Vec3d& east, Vec3d& north, Vec3d& up) noexcept;

// ─────────────────────────────────────────────────────────────────────────────
// Build a full SensorInfo struct from geodetic coordinates
// ─────────────────────────────────────────────────────────────────────────────
SensorInfo make_sensor_info(const std::string& id,
                            double lat_deg, double lon_deg, double alt_m);

// ─────────────────────────────────────────────────────────────────────────────
// Rotation matrix from ECEF to local ENU (3×3)
// ─────────────────────────────────────────────────────────────────────────────
Mat3d ecef_to_enu_rotation(double lat_deg, double lon_deg) noexcept;

}  // namespace geodesy
}  // namespace mtt
