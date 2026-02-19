#pragma once
#include "types.hpp"
#include "track.hpp"
#include "object_pool.hpp"
#include "ukf.hpp"
#include "gating.hpp"
#include "assignment.hpp"
#include "ws3d.hpp"
#include "sensor_manager.hpp"
#include <vector>
#include <unordered_map>
#include <mutex>
#include <algorithm>
#include <shared_mutex>


namespace mtt {

// Forward-declare FullRawMeasurement (defined in input_parser.hpp)
struct FullRawMeasurement;

// ─────────────────────────────────────────────────────────────────────────────
// Pending bearing rays, buffered per timestamp bucket for triangulation.
// ─────────────────────────────────────────────────────────────────────────────
struct PendingRays {
    double                   timestamp = 0.0;
    std::vector<BearingRay>  rays;
    std::vector<SensorID>    sensor_ids;
    std::vector<Measurement> meas;
};

// ─────────────────────────────────────────────────────────────────────────────
// Track-Oriented Multiple Hypothesis Tracker
// ─────────────────────────────────────────────────────────────────────────────
class TOMHT {
public:
    static constexpr double SYNC_WINDOW = 0.05;   // 50 ms triangulation window

    explicit TOMHT(MHTConfig cfg, SensorManager& sm)
        : cfg_(std::move(cfg))
        , sm_(sm)
        , ukf_(cfg_)
        , gate_(cfg_.gate_chi2_threshold) {
        ukf_.psd = cfg_.psd_default;
    }

    // ─── Primary ingest entry point ───────────────────────────────────────────
    void ingest(const Measurement& meas, const SensorInfo& sensor) {
        buffer_and_maybe_birth(meas, sensor);
        process_measurement(meas, sensor);
        ++ingest_count_;
        if (ingest_count_ % 200 == 0) prune_tracks();
    }

    // ─── Thread-safe output snapshot ──────────────────────────────────────────
    std::vector<Track::Snapshot> get_outputs() const {
        std::shared_lock lk(tracks_mx_);
        std::vector<Track::Snapshot> out;
        out.reserve(active_tracks_.size());
        for (const Track* t : active_tracks_)
            if (t->is_output()) out.push_back(t->snapshot());
        return out;
    }

    std::size_t active_track_count() const {
        std::shared_lock lk(tracks_mx_);
        return active_tracks_.size();
    }

private:
    MHTConfig       cfg_;
    SensorManager&  sm_;
    UKF             ukf_;
    CoarseGate      gate_;

    ObjectPool<Track, 10'000>    pool_;
    mutable std::shared_mutex    tracks_mx_;
    std::vector<Track*>          active_tracks_;

    std::mutex                               ray_mx_;
    std::unordered_map<int64_t, PendingRays> ray_buffer_;

    uint64_t ingest_count_ = 0;

    // ─── Buffer ray and trigger triangulation when enough sensors ─────────────
    bool buffer_and_maybe_birth(const Measurement& meas, const SensorInfo& sensor) {
        const double s2 = cfg_.meas_noise_deg2 * DEG2RAD * DEG2RAD;
        BearingRay r{
            sensor.ecef,
            meas.z(0), meas.z(1),
            s2, s2,
            sensor.east_ecef, sensor.north_ecef, sensor.up_ecef
        };

        const int64_t bucket = static_cast<int64_t>(meas.timestamp / SYNC_WINDOW);
        std::lock_guard lk(ray_mx_);
        auto& pb = ray_buffer_[bucket];
        pb.timestamp = meas.timestamp;

        // One ray per sensor per window
        for (const auto& sid : pb.sensor_ids)
            if (sid == sensor.id) return false;

        pb.rays.push_back(r);
        pb.sensor_ids.push_back(sensor.id);
        pb.meas.push_back(meas);

        if (static_cast<int>(pb.rays.size()) >= cfg_.min_sensors_for_init) {
            attempt_track_birth(pb);
            ray_buffer_.erase(bucket);
            return true;
        }
        return false;
    }

    // ─── WS3D birth ───────────────────────────────────────────────────────────
    void attempt_track_birth(const PendingRays& pb) {
        WS3DResult ws = WS3DTriangulator::triangulate(pb.rays);
        if (!ws.valid) return;

        {
            std::shared_lock lk(tracks_mx_);
            for (const Track* t : active_tracks_)
                if ((t->ukf.x.head<3>() - ws.position_ecef).norm() < 300.0)
                    return;
        }

        Track* t = pool_.acquire();
        if (!t) return;
        t->initialise(ws.position_ecef, ws.cov_pos, pb.timestamp, cfg_);

        std::unique_lock lk(tracks_mx_);
        active_tracks_.push_back(t);
    }

    // ─── Gating → Assignment → UKF update ─────────────────────────────────────
    void process_measurement(const Measurement& meas, const SensorInfo& sensor) {
        std::vector<Track*> tracks;
        {
            std::shared_lock lk(tracks_mx_);
            tracks = active_tracks_;
        }
        if (tracks.empty()) return;

        auto gated = gate_.filter(meas, sensor, tracks, ukf_);
        if (gated.empty()) return;

        const int n = static_cast<int>(gated.size());
        std::vector<std::vector<double>> cost(1, std::vector<double>(n));
        for (int j = 0; j < n; ++j) cost[0][j] = gated[j].maha2;

        auto asgn = HungarianSolver::solve(cost, cfg_.gate_chi2_threshold * 1.5);
        if (!asgn.empty() && asgn[0] >= 0 && asgn[0] < n) {
            Track* best = gated[asgn[0]].track;
            best->apply_update(meas, sensor, ukf_);
            best->update_state(cfg_);
        }
    }

    // ─── Prune deleted tracks ─────────────────────────────────────────────────
    void prune_tracks() {
        std::unique_lock lk(tracks_mx_);
        std::vector<Track*> survivors;
        survivors.reserve(active_tracks_.size());
        for (Track* t : active_tracks_) {
            if (t->is_alive()) survivors.push_back(t);
            else pool_.release(t);
        }
        active_tracks_ = std::move(survivors);
    }
};

}  // namespace mtt
