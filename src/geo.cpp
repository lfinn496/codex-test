#include "mtt/geo.hpp"

#include <cmath>

namespace {
constexpr double kWgs84A = 6378137.0;
constexpr double kWgs84F = 1.0 / 298.257223563;
constexpr double kPi = 3.14159265358979323846;

double deg2rad(double deg) {
    return deg * (kPi / 180.0);
}
} // namespace

namespace mtt {

EcefPoint wgs84_to_ecef(const SensorPoseWgs84& lla) {
    const double lat = deg2rad(lla.latitude_deg);
    const double lon = deg2rad(lla.longitude_deg);

    const double e2 = kWgs84F * (2.0 - kWgs84F);
    const double sin_lat = std::sin(lat);
    const double cos_lat = std::cos(lat);
    const double sin_lon = std::sin(lon);
    const double cos_lon = std::cos(lon);

    const double N = kWgs84A / std::sqrt(1.0 - e2 * sin_lat * sin_lat);

    return {
        (N + lla.altitude_m) * cos_lat * cos_lon,
        (N + lla.altitude_m) * cos_lat * sin_lon,
        (N * (1.0 - e2) + lla.altitude_m) * sin_lat,
    };
}

std::array<double, 3> los_enu_from_az_el_deg(double azimuth_deg, double elevation_deg) {
    const double az = deg2rad(azimuth_deg);
    const double el = deg2rad(elevation_deg);

    const double cos_el = std::cos(el);
    const double east = cos_el * std::sin(az);
    const double north = cos_el * std::cos(az);
    const double up = std::sin(el);
    return {east, north, up};
}

std::array<double, 3> enu_to_ecef_direction(const SensorPoseWgs84& sensor, const std::array<double, 3>& enu_dir) {
    const double lat = deg2rad(sensor.latitude_deg);
    const double lon = deg2rad(sensor.longitude_deg);

    const double sin_lat = std::sin(lat);
    const double cos_lat = std::cos(lat);
    const double sin_lon = std::sin(lon);
    const double cos_lon = std::cos(lon);

    // ENU basis in ECEF coordinates.
    const std::array<double, 3> e{-sin_lon, cos_lon, 0.0};
    const std::array<double, 3> n{-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat};
    const std::array<double, 3> u{cos_lat * cos_lon, cos_lat * sin_lon, sin_lat};

    return {
        enu_dir[0] * e[0] + enu_dir[1] * n[0] + enu_dir[2] * u[0],
        enu_dir[0] * e[1] + enu_dir[1] * n[1] + enu_dir[2] * u[1],
        enu_dir[0] * e[2] + enu_dir[1] * n[2] + enu_dir[2] * u[2],
    };
}

} // namespace mtt
