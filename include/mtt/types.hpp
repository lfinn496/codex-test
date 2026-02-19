#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace mtt {

struct SensorPoseWgs84 {
    double latitude_deg{0.0};
    double longitude_deg{0.0};
    double altitude_m{0.0};
};

struct Measurement {
    double timestamp_s{0.0};
    double azimuth_deg{0.0};
    double elevation_deg{0.0};
    SensorPoseWgs84 sensor;
};

struct Covariance2x2 {
    // Units: deg^2. Symmetric matrix.
    double xx{1.0};
    double xy{0.0};
    double yy{1.0};
};

struct TrackerConfig {
    Covariance2x2 measurement_cov_deg2{};
    double process_noise_accel_mps2{4.0};
    double survival_probability{0.997};
    double detection_probability{0.9};
    double clutter_intensity{1e-3};
    double gating_threshold_mahalanobis{9.21}; // chi^2(2 dof, 0.99)
    std::size_t max_hypotheses{2000};
    std::size_t max_tracks{5000};
};

struct TrackState {
    // Cartesian ECEF state [x y z vx vy vz].
    std::array<double, 6> x{};
    std::array<double, 36> P{};
};

struct Track {
    std::uint64_t id{0};
    double existence_probability{0.0};
    TrackState state{};
    std::uint32_t age{0};
    std::uint32_t miss_streak{0};
};

struct UpdateBatch {
    double timestamp_s{0.0};
    std::vector<Measurement> measurements;
};

} // namespace mtt
