#include "mtt/geodesy.hpp"
#include <cmath>

namespace mtt {
namespace geodesy {

// ─────────────────────────────────────────────────────────────────────────────
Vec3d lla_to_ecef(double lat_deg, double lon_deg, double alt_m) noexcept {
    const double lat = deg2rad(lat_deg);
    const double lon = deg2rad(lon_deg);
    const double N   = prime_vertical_radius(lat);
    const double cos_lat = std::cos(lat);
    const double sin_lat = std::sin(lat);
    const double cos_lon = std::cos(lon);
    const double sin_lon = std::sin(lon);
    Vec3d result;
    result << (N + alt_m) * cos_lat * cos_lon,
              (N + alt_m) * cos_lat * sin_lon,
              (N * (1.0 - wgs84::e2) + alt_m) * sin_lat;
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
void ecef_to_lla(const Vec3d& ecef,
                 double& lat_deg, double& lon_deg, double& alt_m) noexcept {
    const double x = ecef.x(), y = ecef.y(), z = ecef.z();
    const double p   = std::sqrt(x*x + y*y);
    lon_deg = rad2deg(std::atan2(y, x));

    // Bowring's iterative method
    double theta = std::atan2(z * wgs84::a, p * wgs84::b);
    double lat   = std::atan2(z + wgs84::ep2 * wgs84::b * std::pow(std::sin(theta),3),
                              p - wgs84::e2  * wgs84::a * std::pow(std::cos(theta),3));
    for (int i = 0; i < 5; ++i) {
        const double N_i = prime_vertical_radius(lat);
        const double lat_new = std::atan2(z + wgs84::e2 * N_i * std::sin(lat), p);
        if (std::abs(lat_new - lat) < 1e-12) { lat = lat_new; break; }
        lat = lat_new;
    }
    const double N = prime_vertical_radius(lat);
    const double cos_lat = std::cos(lat);
    alt_m   = (std::abs(cos_lat) > 1e-9) ? (p / cos_lat - N)
                                          : (std::abs(z) / std::sin(lat) - N*(1.0-wgs84::e2));
    lat_deg = rad2deg(lat);
}

// ─────────────────────────────────────────────────────────────────────────────
void compute_enu_axes(double lat_deg, double lon_deg,
                      Vec3d& east, Vec3d& north, Vec3d& up) noexcept {
    const double lat = deg2rad(lat_deg);
    const double lon = deg2rad(lon_deg);
    const double cl = std::cos(lat), sl = std::sin(lat);
    const double co = std::cos(lon), so = std::sin(lon);
    // East  = [-sin(lon),  cos(lon),      0       ]
    east << -so,   co,   0.0 ;
    // North = [-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)]
    north << -sl*co, -sl*so,  cl  ;
    // Up    = [ cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat)]
    up << cl*co,   cl*so,  sl  ;
}

// ─────────────────────────────────────────────────────────────────────────────
SensorInfo make_sensor_info(const std::string& id,
                            double lat_deg, double lon_deg, double alt_m) {
    SensorInfo s;
    s.id      = id;
    s.lat_deg = lat_deg;
    s.lon_deg = lon_deg;
    s.alt_m   = alt_m;
    s.ecef    = lla_to_ecef(lat_deg, lon_deg, alt_m);
    compute_enu_axes(lat_deg, lon_deg, s.east_ecef, s.north_ecef, s.up_ecef);
    return s;
}

// ─────────────────────────────────────────────────────────────────────────────
Mat3d ecef_to_enu_rotation(double lat_deg, double lon_deg) noexcept {
    Vec3d e, n, u;
    compute_enu_axes(lat_deg, lon_deg, e, n, u);
    Mat3d R;
    R(0,0)=e(0); R(0,1)=e(1); R(0,2)=e(2);   // East  row
    R(1,0)=n(0); R(1,1)=n(1); R(1,2)=n(2);   // North row
    R(2,0)=u(0); R(2,1)=u(1); R(2,2)=u(2);   // Up    row
    return R;
}

}  // namespace geodesy
}  // namespace mtt
