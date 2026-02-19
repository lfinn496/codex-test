#include "mtt/pmbm_tracker.hpp"
#include "mtt/geo.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

double dot3(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<double, 3> sub3(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

std::array<double, 3> add3(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

std::array<double, 3> scale3(const std::array<double, 3>& a, double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}

double norm3(const std::array<double, 3>& v) {
    return std::sqrt(dot3(v, v));
}

double wrap_deg(double angle_deg) {
    double a = std::fmod(angle_deg + 180.0, 360.0);
    if (a < 0.0) {
        a += 360.0;
    }
    return a - 180.0;
}

double rad2deg(double rad) {
    return rad * (180.0 / 3.14159265358979323846);
}

std::array<double, 6> constant_velocity_predict(const std::array<double, 6>& x, double dt_s) {
    return {
        x[0] + dt_s * x[3],
        x[1] + dt_s * x[4],
        x[2] + dt_s * x[5],
        x[3],
        x[4],
        x[5],
    };
}

struct SensorRay {
    std::size_t meas_index{0};
    std::array<double, 3> origin{};
    std::array<double, 3> direction{};
    mtt::SensorPoseWgs84 sensor;
};

struct CandidatePoint {
    std::array<double, 3> p{};
    double ray_separation_m{0.0};
};

bool triangulate_two_rays(const SensorRay& a, const SensorRay& b, CandidatePoint& out) {
    const auto w0 = sub3(a.origin, b.origin);
    const double aa = dot3(a.direction, a.direction);
    const double bb = dot3(a.direction, b.direction);
    const double cc = dot3(b.direction, b.direction);
    const double dd = dot3(a.direction, w0);
    const double ee = dot3(b.direction, w0);

    const double denom = aa * cc - bb * bb;
    if (std::abs(denom) < 1e-9) {
        return false;
    }

    const double ta = (bb * ee - cc * dd) / denom;
    const double tb = (aa * ee - bb * dd) / denom;

    const auto pa = add3(a.origin, scale3(a.direction, ta));
    const auto pb = add3(b.origin, scale3(b.direction, tb));

    const auto mid = scale3(add3(pa, pb), 0.5);
    out.p = mid;
    out.ray_separation_m = norm3(sub3(pa, pb));
    return true;
}

std::pair<double, double> predict_az_el_deg(const mtt::SensorPoseWgs84& sensor, const std::array<double, 3>& target_ecef) {
    constexpr double kPi = 3.14159265358979323846;
    const double lat = sensor.latitude_deg * (kPi / 180.0);
    const double lon = sensor.longitude_deg * (kPi / 180.0);

    const double sin_lat = std::sin(lat);
    const double cos_lat = std::cos(lat);
    const double sin_lon = std::sin(lon);
    const double cos_lon = std::cos(lon);

    const auto sensor_ecef = mtt::wgs84_to_ecef(sensor);
    const std::array<double, 3> r{target_ecef[0] - sensor_ecef.x_m, target_ecef[1] - sensor_ecef.y_m, target_ecef[2] - sensor_ecef.z_m};

    const double east = -sin_lon * r[0] + cos_lon * r[1];
    const double north = -sin_lat * cos_lon * r[0] - sin_lat * sin_lon * r[1] + cos_lat * r[2];
    const double up = cos_lat * cos_lon * r[0] + cos_lat * sin_lon * r[1] + sin_lat * r[2];

    const double az = rad2deg(std::atan2(east, north));
    const double el = rad2deg(std::atan2(up, std::hypot(east, north)));
    return {az, el};
}

} // namespace

namespace mtt {

struct PmbmTracker::Impl {
    explicit Impl(TrackerConfig cfg)
        : config(std::move(cfg)) {}

    TrackerConfig config;
    double current_time_s{0.0};
    bool initialized{false};
    std::uint64_t next_track_id{1};
    std::vector<Track> tracks;
    double last_update_time_s{0.0};
    bool has_update_time{false};

    static constexpr double kMaxPairSeparationM = 250.0;
    static constexpr double kResidualAzElDeg = 1.0;
    static constexpr double kTrackAssociationGateM = 350.0;
    static constexpr double kClusterDistanceM = 150.0;

    void predict(double timestamp_s) {
        if (!initialized) {
            current_time_s = timestamp_s;
            initialized = true;
            return;
        }

        const double dt = std::max(0.0, timestamp_s - current_time_s);
        current_time_s = timestamp_s;

        for (auto& trk : tracks) {
            trk.state.x = constant_velocity_predict(trk.state.x, dt);
            trk.existence_probability *= config.survival_probability;
            trk.age += 1;
            trk.miss_streak += 1;
        }
    }

    std::vector<SensorRay> build_rays(const std::vector<Measurement>& measurements) const {
        std::vector<SensorRay> rays;
        rays.reserve(measurements.size());
        for (std::size_t i = 0; i < measurements.size(); ++i) {
            const auto& z = measurements[i];
            const auto sensor_ecef = wgs84_to_ecef(z.sensor);
            const auto enu = los_enu_from_az_el_deg(z.azimuth_deg, z.elevation_deg);
            const auto ecef_dir = enu_to_ecef_direction(z.sensor, enu);
            rays.push_back({i, {sensor_ecef.x_m, sensor_ecef.y_m, sensor_ecef.z_m}, ecef_dir, z.sensor});
        }
        return rays;
    }

    std::vector<std::array<double, 3>> triangulate_candidates(const std::vector<Measurement>& measurements) const {
        const auto rays = build_rays(measurements);
        std::vector<std::array<double, 3>> candidates;

        for (std::size_t i = 0; i < rays.size(); ++i) {
            for (std::size_t j = i + 1; j < rays.size(); ++j) {
                const auto& a = rays[i];
                const auto& b = rays[j];
                if (a.sensor.latitude_deg == b.sensor.latitude_deg && a.sensor.longitude_deg == b.sensor.longitude_deg &&
                    a.sensor.altitude_m == b.sensor.altitude_m) {
                    continue;
                }

                CandidatePoint cand;
                if (!triangulate_two_rays(a, b, cand)) {
                    continue;
                }
                if (cand.ray_separation_m > kMaxPairSeparationM) {
                    continue;
                }

                const auto [az_a, el_a] = predict_az_el_deg(a.sensor, cand.p);
                const auto [az_b, el_b] = predict_az_el_deg(b.sensor, cand.p);
                const auto& za = measurements[a.meas_index];
                const auto& zb = measurements[b.meas_index];

                if (std::abs(wrap_deg(za.azimuth_deg - az_a)) > kResidualAzElDeg || std::abs(za.elevation_deg - el_a) > kResidualAzElDeg ||
                    std::abs(wrap_deg(zb.azimuth_deg - az_b)) > kResidualAzElDeg || std::abs(zb.elevation_deg - el_b) > kResidualAzElDeg) {
                    continue;
                }

                candidates.push_back(cand.p);
            }
        }

        return cluster_candidates(candidates);
    }

    static std::vector<std::array<double, 3>> cluster_candidates(const std::vector<std::array<double, 3>>& candidates) {
        std::vector<std::array<double, 3>> centroids;
        std::vector<std::size_t> counts;

        for (const auto& c : candidates) {
            bool assigned = false;
            for (std::size_t i = 0; i < centroids.size(); ++i) {
                if (norm3(sub3(c, centroids[i])) < kClusterDistanceM) {
                    const double n = static_cast<double>(counts[i]);
                    centroids[i] = scale3(add3(scale3(centroids[i], n), c), 1.0 / (n + 1.0));
                    counts[i] += 1;
                    assigned = true;
                    break;
                }
            }
            if (!assigned) {
                centroids.push_back(c);
                counts.push_back(1);
            }
        }
        return centroids;
    }

    void update_track_from_point(Track& trk, const std::array<double, 3>& p_meas, double dt_s) {
        std::array<double, 3> prev_pos{trk.state.x[0], trk.state.x[1], trk.state.x[2]};
        std::array<double, 3> prev_vel{trk.state.x[3], trk.state.x[4], trk.state.x[5]};

        const double alpha = 0.65;
        const double beta = 0.4;

        const auto pos_new = add3(scale3(prev_pos, 1.0 - alpha), scale3(p_meas, alpha));
        std::array<double, 3> vel_new = prev_vel;
        if (dt_s > 1e-3) {
            const auto fd_vel = scale3(sub3(pos_new, prev_pos), 1.0 / dt_s);
            vel_new = add3(scale3(prev_vel, 1.0 - beta), scale3(fd_vel, beta));
        }

        trk.state.x = {pos_new[0], pos_new[1], pos_new[2], vel_new[0], vel_new[1], vel_new[2]};
        trk.existence_probability = std::min(0.999, trk.existence_probability + 0.2);
        trk.miss_streak = 0;
    }

    void spawn_birth_track(const std::array<double, 3>& p) {
        Track trk;
        trk.id = next_track_id++;
        trk.existence_probability = 0.55;
        trk.state.x = {p[0], p[1], p[2], 0.0, 0.0, 0.0};
        trk.age = 1;
        tracks.push_back(trk);
    }

    void update(const UpdateBatch& batch) {
        const auto candidates = triangulate_candidates(batch.measurements);
        const double dt = has_update_time ? std::max(0.0, batch.timestamp_s - last_update_time_s) : 0.0;

        std::vector<bool> used_track(tracks.size(), false);
        for (const auto& c : candidates) {
            double best_dist = std::numeric_limits<double>::max();
            std::size_t best_idx = 0;
            bool found = false;
            for (std::size_t i = 0; i < tracks.size(); ++i) {
                if (used_track[i]) {
                    continue;
                }
                const std::array<double, 3> p_pred{tracks[i].state.x[0], tracks[i].state.x[1], tracks[i].state.x[2]};
                const double d = norm3(sub3(c, p_pred));
                if (d < kTrackAssociationGateM && d < best_dist) {
                    best_dist = d;
                    best_idx = i;
                    found = true;
                }
            }

            if (found) {
                update_track_from_point(tracks[best_idx], c, dt);
                used_track[best_idx] = true;
            } else {
                spawn_birth_track(c);
            }
        }

        for (auto& trk : tracks) {
            if (trk.miss_streak > 0) {
                trk.existence_probability *= 0.92;
            }
        }

        tracks.erase(std::remove_if(tracks.begin(), tracks.end(), [](const Track& trk) {
                         return trk.existence_probability < 0.1 || trk.miss_streak > 15;
                     }),
            tracks.end());

        if (tracks.size() > config.max_tracks) {
            tracks.erase(tracks.begin(), tracks.begin() + static_cast<std::ptrdiff_t>(tracks.size() - config.max_tracks));
        }

        last_update_time_s = batch.timestamp_s;
        has_update_time = true;
    }
};

PmbmTracker::PmbmTracker(TrackerConfig config)
    : impl_(std::make_unique<Impl>(std::move(config))) {}

PmbmTracker::~PmbmTracker() = default;

void PmbmTracker::predict_to(double timestamp_s) {
    impl_->predict(timestamp_s);
}

void PmbmTracker::update(const UpdateBatch& batch) {
    impl_->update(batch);
}

std::vector<Track> PmbmTracker::tracks() const {
    return impl_->tracks;
}

} // namespace mtt
