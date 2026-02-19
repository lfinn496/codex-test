#pragma once
#include "types.hpp"
#include "blocking_queue.hpp"
#include "sensor_manager.hpp"
#include <istream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// Input format (space or comma separated):
//   timestamp  azimuth_deg  elevation_deg  sensor_lat_deg  sensor_lon_deg  sensor_alt_m
//
// Lines starting with '#' are treated as comments and ignored.
// ─────────────────────────────────────────────────────────────────────────────

// Extend RawMeasurement to carry sensor geodetic coords for registration
struct FullRawMeasurement : public RawMeasurement {
    double sensor_lat_deg;
    double sensor_lon_deg;
    double sensor_alt_m;
};

class InputParser {
public:
    explicit InputParser(SensorManager& sm) : sm_(sm) {}

    // Parse one line; returns true if parsing succeeded.
    bool parse_line(const std::string& line, FullRawMeasurement& out) const {
        if (line.empty() || line[0] == '#') return false;

        // Replace commas with spaces for uniform tokenisation
        std::string normalized = line;
        for (char& c : normalized) if (c == ',') c = ' ';

        std::istringstream iss(normalized);
        double ts, az, el, lat, lon, alt;
        if (!(iss >> ts >> az >> el >> lat >> lon >> alt)) return false;

        out.timestamp     = ts;
        out.az_deg        = az;
        out.el_deg        = el;
        out.sensor_lat_deg= lat;
        out.sensor_lon_deg= lon;
        out.sensor_alt_m  = alt;
        out.sensor_id     = SensorManager::make_id(lat, lon, alt);
        return true;
    }

    // Stream all lines from istream into the queue.
    // Blocks until EOF or until the queue signals stop.
    void stream_into(std::istream& in,
                     BlockingQueue<FullRawMeasurement>& q,
                     bool verbose = false) {
        std::string line;
        uint64_t line_no = 0;
        uint64_t parsed  = 0;
        while (std::getline(in, line)) {
            ++line_no;
            FullRawMeasurement m;
            if (!parse_line(line, m)) continue;

            // Register sensor (cached after first encounter)
            sm_.get_or_register(m.sensor_lat_deg, m.sensor_lon_deg, m.sensor_alt_m);

            q.push(m);
            ++parsed;
            if (verbose && parsed % 10'000 == 0)
                std::cerr << "[InputParser] Parsed " << parsed
                          << " measurements (line " << line_no << ")\n";
            if (q.is_stopped()) break;
        }
        q.stop();
        if (verbose)
            std::cerr << "[InputParser] Done. Total lines=" << line_no
                      << " measurements=" << parsed << "\n";
    }

private:
    SensorManager& sm_;
};

}  // namespace mtt
