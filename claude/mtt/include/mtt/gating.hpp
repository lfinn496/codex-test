#pragma once
#include "types.hpp"
#include "track.hpp"
#include "ukf.hpp"
#include "sensor_manager.hpp"
#include <vector>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// Coarse Gating (Tier 1 of the 3-tier data association pipeline).
//
// For each incoming measurement, computes the Mahalanobis² distance to every
// active track's predicted measurement. Returns only the (track, maha²) pairs
// that fall within the chi-squared gate.
//
// Chi-squared thresholds for 2 DoF (azimuth + elevation):
//   90%  → 4.605
//   95%  → 5.991
//   99%  → 9.210
//   99.9%→ 13.816
// ─────────────────────────────────────────────────────────────────────────────
struct GatedPair {
    Track*  track;
    double  maha2;
};

class CoarseGate {
public:
    explicit CoarseGate(double chi2_threshold)
        : threshold_(chi2_threshold) {}

    // For one measurement, return all tracks that pass the gate.
    std::vector<GatedPair> filter(const Measurement& meas,
                                  const SensorInfo&  sensor,
                                  const std::vector<Track*>& tracks,
                                  const UKF& filter) const {
        std::vector<GatedPair> passed;
        passed.reserve(tracks.size());
        for (Track* t : tracks) {
            const double maha2 = filter.gating_distance(t->ukf, meas, sensor);
            if (maha2 <= threshold_) {
                passed.push_back({t, maha2});
            }
        }
        return passed;
    }

    double threshold() const noexcept { return threshold_; }
    void   set_threshold(double v) noexcept { threshold_ = v; }

private:
    double threshold_;
};

}  // namespace mtt
