#pragma once
#include "types.hpp"
#include "geodesy.hpp"
#include <unordered_map>
#include <mutex>
#include <sstream>
#include <iomanip>
#include <shared_mutex>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// Thread-safe sensor registry.
// The first time a (lat, lon, alt) triple is encountered in the input stream,
// the ECEF position and ENU axes are computed and cached permanently.
// Subsequent lookups cost only a hash-table read.
// ─────────────────────────────────────────────────────────────────────────────
class SensorManager {
public:
    // Tolerance for considering two geodetic coordinates the same sensor.
    static constexpr double LAT_LON_TOL_DEG = 1e-6;   // ~0.1 m
    static constexpr double ALT_TOL_M       = 0.1;

    // Build a stable string key from geodetic coords.
    static SensorID make_id(double lat_deg, double lon_deg, double alt_m) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(6)
           << lat_deg << "," << lon_deg << "," << std::setprecision(1) << alt_m;
        return ss.str();
    }

    // Register (if not already present) and return the SensorInfo reference.
    const SensorInfo& get_or_register(double lat_deg, double lon_deg, double alt_m) {
        const SensorID id = make_id(lat_deg, lon_deg, alt_m);
        {
            std::shared_lock lk(mx_);
            auto it = sensors_.find(id);
            if (it != sensors_.end()) return it->second;
        }
        // Not found – compute and insert (write lock).
        SensorInfo info = geodesy::make_sensor_info(id, lat_deg, lon_deg, alt_m);
        std::unique_lock lk(mx_);
        return sensors_.emplace(id, std::move(info)).first->second;
    }

    const SensorInfo* find(const SensorID& id) const {
        std::shared_lock lk(mx_);
        auto it = sensors_.find(id);
        if (it == sensors_.end()) return nullptr;
        return &it->second;
    }

    std::size_t size() const {
        std::shared_lock lk(mx_);
        return sensors_.size();
    }

private:
    mutable std::shared_mutex mx_;
    std::unordered_map<SensorID, SensorInfo> sensors_;
};

}  // namespace mtt
