#include "mtt/pmbm_tracker.hpp"
#include "mtt/geo.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>

namespace {

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
        }
    }

    void spawn_birth_track(const Measurement& z) {
        const auto sensor_ecef = wgs84_to_ecef(z.sensor);
        const auto dir_enu = los_enu_from_az_el_deg(z.azimuth_deg, z.elevation_deg);
        const auto dir_ecef = enu_to_ecef_direction(z.sensor, dir_enu);

        // Simple birth prior along line-of-sight, seeded at nominal range.
        constexpr double nominal_range_m = 2000.0;

        Track trk;
        trk.id = next_track_id++;
        trk.existence_probability = 0.2;
        trk.state.x = {
            sensor_ecef.x_m + nominal_range_m * dir_ecef[0],
            sensor_ecef.y_m + nominal_range_m * dir_ecef[1],
            sensor_ecef.z_m + nominal_range_m * dir_ecef[2],
            0.0,
            0.0,
            0.0,
        };
        trk.age = 1;
        tracks.push_back(trk);
    }

    void lightweight_update(const UpdateBatch& batch) {
        // Placeholder update strategy: spawn a Bernoulli birth for each
        // measurement. Production system should replace this with full PMBM
        // update, Murty k-best assignment, and hypothesis pruning.
        for (const auto& z : batch.measurements) {
            spawn_birth_track(z);
        }

        if (tracks.size() > config.max_tracks) {
            const auto to_remove = tracks.size() - config.max_tracks;
            tracks.erase(tracks.begin(), tracks.begin() + static_cast<std::ptrdiff_t>(to_remove));
        }

        for (auto& trk : tracks) {
            trk.existence_probability = std::min(0.999, trk.existence_probability + 0.05);
        }

        tracks.erase(std::remove_if(tracks.begin(), tracks.end(), [](const Track& trk) {
                         return trk.existence_probability < 1e-3;
                     }),
            tracks.end());
    }
};

PmbmTracker::PmbmTracker(TrackerConfig config)
    : impl_(std::make_unique<Impl>(std::move(config))) {}

PmbmTracker::~PmbmTracker() = default;

void PmbmTracker::predict_to(double timestamp_s) {
    impl_->predict(timestamp_s);
}

void PmbmTracker::update(const UpdateBatch& batch) {
    impl_->lightweight_update(batch);
}

std::vector<Track> PmbmTracker::tracks() const {
    return impl_->tracks;
}

} // namespace mtt
