#pragma once
#include "types.hpp"
#include "ukf.hpp"
#include <atomic>
#include <deque>
#include <mutex>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// Single-track hypothesis node.
//
// Each Track owns:
//  - a UKF state (mean + covariance)
//  - a log-likelihood score for TOMHT pruning
//  - hit/miss counters for lifecycle management
//  - a fixed-depth innovation history (for N-scan pruning)
// ─────────────────────────────────────────────────────────────────────────────
class Track {
public:
    static std::atomic<TrackID> id_counter;

    TrackID    id          = INVALID_TRACK_ID;
    TrackState state       = TrackState::Tentative;

    UKFState   ukf;

    // Log-likelihood score (TOMHT)
    double     score       = 0.0;      // log P(track | observations)

    // Consecutive hit/miss counters
    int        hits        = 0;
    int        misses      = 0;
    int        total_hits  = 0;

    // Innovation history (Mahalanobis² values, size-limited to n_scan_depth)
    std::deque<double> innov_history;

    // Mutex protecting this track from concurrent UKF updates
    mutable std::mutex mtx;

    // ─── Reset (called by ObjectPool on release) ──────────────────────────────
    void reset() noexcept {
        id     = INVALID_TRACK_ID;
        state  = TrackState::Tentative;
        ukf.x  = Vec6d::Zero();
        ukf.P  = Mat6d::Identity();
        ukf.timestamp = 0.0;
        score  = 0.0;
        hits   = 0;
        misses = 0;
        total_hits = 0;
        innov_history.clear();
    }

    // ─── Initialise from WS3D result ──────────────────────────────────────────
    void initialise(const Vec3d& pos_ecef, const Mat3d& cov_pos,
                    double timestamp, const MHTConfig& cfg) {
        id = ++id_counter;
        state = TrackState::Tentative;

        ukf.x = Vec6d::Zero();
        ukf.x.head<3>() = pos_ecef;
        // Velocity initialised at zero; covariance inflated

        ukf.P = Mat6d::Zero();
        ukf.P.topLeftCorner<3,3>()     = cov_pos;
        ukf.P.bottomRightCorner<3,3>() = Mat3d::Identity() * cfg.init_vel_cov;
        ukf.timestamp = timestamp;

        score  = 0.0;
        hits   = 0;
        misses = 0;
        total_hits = 0;
    }

    // ─── Apply UKF predict+update (sequential asynchronous update) ─────────────
    double apply_update(const Measurement& meas, const SensorInfo& sensor,
                        const UKF& filter) {
        std::lock_guard lk(mtx);

        // Prediction step: propagate to measurement timestamp
        double dt = meas.timestamp - ukf.timestamp;
        if (dt < 0.0) dt = 0.0;
        if (dt > 0.0) filter.predict(ukf, dt);
        ukf.timestamp = meas.timestamp;

        // Measurement update
        const double maha2 = filter.update(ukf, meas, sensor);
        ukf.timestamp = meas.timestamp;

        // Score update: log-likelihood contribution
        // Approximated as -0.5 * maha2 for chi2 distribution
        score += -0.5 * maha2;

        // Hit/miss bookkeeping
        ++hits;
        ++total_hits;
        misses = 0;
        innov_history.push_back(maha2);
        return maha2;
    }

    // ─── Record a missed update ───────────────────────────────────────────────
    void apply_miss(double timestamp, const UKF& filter) {
        std::lock_guard lk(mtx);
        double dt = timestamp - ukf.timestamp;
        if (dt > 0.0) filter.predict(ukf, dt);
        ukf.timestamp = timestamp;
        ++misses;
        // Score penalty for missed detection
        score -= 1.5;
        innov_history.push_back(0.0);
    }

    // ─── Lifecycle decisions ──────────────────────────────────────────────────
    void update_state(const MHTConfig& cfg) noexcept {
        // Trim innovation history
        while (static_cast<int>(innov_history.size()) > cfg.n_scan_depth)
            innov_history.pop_front();

        if (state == TrackState::Deleted) return;

        if (state == TrackState::Tentative) {
            if (hits >= cfg.confirm_hits)
                state = TrackState::Confirmed;
        }
        if (misses >= cfg.deletion_misses || score <= cfg.deletion_score_floor)
            state = TrackState::Deleted;
        if (misses > 0 && state == TrackState::Confirmed)
            state = TrackState::Coasting;
        if (misses == 0 && state == TrackState::Coasting)
            state = TrackState::Confirmed;
    }

    bool is_alive() const noexcept {
        return state != TrackState::Deleted;
    }

    bool is_output() const noexcept {
        return state == TrackState::Confirmed || state == TrackState::Coasting;
    }

    // ─── Thread-safe state snapshot for output ────────────────────────────────
    struct Snapshot {
        TrackID    id;
        TrackState state;
        Vec6d      state_vec;
        double     timestamp;
        double     score;
    };

    Snapshot snapshot() const {
        std::lock_guard lk(mtx);
        return { id, state, ukf.x, ukf.timestamp, score };
    }
};

// Global atomic ID counter
inline std::atomic<TrackID> Track::id_counter{0};

}  // namespace mtt
