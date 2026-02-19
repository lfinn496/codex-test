#pragma once

#include "mtt/types.hpp"

#include <array>

namespace mtt {

struct EcefPoint {
    double x_m{0.0};
    double y_m{0.0};
    double z_m{0.0};
};

EcefPoint wgs84_to_ecef(const SensorPoseWgs84& lla);
std::array<double, 3> los_enu_from_az_el_deg(double azimuth_deg, double elevation_deg);
std::array<double, 3> enu_to_ecef_direction(const SensorPoseWgs84& sensor, const std::array<double, 3>& enu_dir);

} // namespace mtt
