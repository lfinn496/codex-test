#include "mtt/pmbm_tracker.hpp"

#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<mtt::Measurement> load_measurements(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open measurements file: " + path);
    }

    std::vector<mtt::Measurement> measurements;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (!std::isdigit(static_cast<unsigned char>(line[0])) && line[0] != '-') {
            continue; // Skip header/comments.
        }

        std::stringstream ss(line);
        std::string token;
        std::vector<double> fields;
        while (std::getline(ss, token, ',')) {
            fields.push_back(std::stod(token));
        }

        if (fields.size() != 6) {
            continue;
        }

        measurements.push_back({
            fields[0],
            fields[1],
            fields[2],
            {fields[3], fields[4], fields[5]},
        });
    }

    return measurements;
}

std::vector<mtt::UpdateBatch> bucket_by_timestamp(const std::vector<mtt::Measurement>& measurements) {
    std::vector<mtt::UpdateBatch> batches;
    constexpr double epsilon = 1e-3;

    for (const auto& z : measurements) {
        if (batches.empty() || std::abs(z.timestamp_s - batches.back().timestamp_s) > epsilon) {
            batches.push_back({z.timestamp_s, {z}});
        } else {
            batches.back().measurements.push_back(z);
        }
    }
    return batches;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: mtt_cli <measurements.csv>\n";
        return 1;
    }

    const auto measurements = load_measurements(argv[1]);
    const auto batches = bucket_by_timestamp(measurements);

    mtt::TrackerConfig cfg;
    cfg.measurement_cov_deg2 = {0.04, 0.0, 0.04};
    mtt::PmbmTracker tracker(cfg);

    for (const auto& batch : batches) {
        tracker.predict_to(batch.timestamp_s);
        tracker.update(batch);
    }

    const auto tracks = tracker.tracks();
    std::cout << "Generated tracks: " << tracks.size() << "\n";
    for (const auto& t : tracks) {
        std::cout << "track_id=" << t.id << " p_exist=" << std::fixed << std::setprecision(3) << t.existence_probability
                  << " ecef_pos_m=[" << t.state.x[0] << ", " << t.state.x[1] << ", " << t.state.x[2] << "]"
                  << " ecef_vel_mps=[" << t.state.x[3] << ", " << t.state.x[4] << ", " << t.state.x[5] << "]\n";
    }

    return 0;
}
