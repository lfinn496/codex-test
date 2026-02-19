#include "mtt/types.hpp"
#include "mtt/geodesy.hpp"
#include "mtt/sensor_manager.hpp"
#include "mtt/observation_model.hpp"
#include "mtt/ukf.hpp"
#include "mtt/ws3d.hpp"
#include "mtt/track.hpp"
#include "mtt/object_pool.hpp"
#include "mtt/gating.hpp"
#include "mtt/assignment.hpp"
#include "mtt/tomht.hpp"
#include "mtt/blocking_queue.hpp"
#include "mtt/input_parser.hpp"
#include "mtt/output_serializer.hpp"

#include <iostream>
#include <fstream>
#include <thread>
#include <atomic>
#include <csignal>
#include <chrono>
#include <string>
#include <stdexcept>

static std::atomic<bool> g_running{true};
static void sig_handler(int) { g_running.store(false, std::memory_order_relaxed); }

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options] [input_file]\n\n"
        << "Input columns (space or comma separated):\n"
        << "  timestamp  azimuth_deg  elevation_deg  sensor_lat  sensor_lon  sensor_alt_m\n\n"
        << "Options:\n"
        << "  --gate N     Chi-sq gate threshold (default 13.816)\n"
        << "  --psd N      Process noise PSD m2/s3 (default 3.0)\n"
        << "  --noise N    Meas noise deg2/axis (default 0.01)\n"
        << "  --confirm N  Hits to confirm track   (default 3)\n"
        << "  --delete N   Misses to delete track  (default 5)\n"
        << "  --nscan N    N-scan pruning depth    (default 5)\n"
        << "  --verbose    Print summaries to stderr\n"
        << "  --out FILE   Output file (default stdout)\n\n"
        << "Output (CSV):\n"
        << "  timestamp,track_id,state,px_m,py_m,pz_m,vx_mps,vy_mps,vz_mps\n";
}

int main(int argc, char** argv) {
    using namespace mtt;
    using namespace std::chrono_literals;

    std::signal(SIGINT,  sig_handler);
    std::signal(SIGTERM, sig_handler);

    MHTConfig cfg;
    std::string input_path, output_path;
    bool verbose = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string {
            if (i+1 >= argc) throw std::invalid_argument("Missing value for " + arg);
            return argv[++i];
        };
        if      (arg == "--gate")    cfg.gate_chi2_threshold = std::stod(next());
        else if (arg == "--psd")     cfg.psd_default         = std::stod(next());
        else if (arg == "--noise")   cfg.meas_noise_deg2     = std::stod(next());
        else if (arg == "--confirm") cfg.confirm_hits        = std::stoi(next());
        else if (arg == "--delete")  cfg.deletion_misses     = std::stoi(next());
        else if (arg == "--nscan")   cfg.n_scan_depth        = std::stoi(next());
        else if (arg == "--verbose") verbose = true;
        else if (arg == "--out")     output_path = next();
        else if (arg == "--help" || arg == "-h") { usage(argv[0]); return 0; }
        else if (!arg.empty() && arg[0] != '-')  input_path = arg;
        else { std::cerr << "Unknown option: " << arg << "\n"; usage(argv[0]); return 1; }
    }

    std::ifstream fin;
    std::istream* in = &std::cin;
    if (!input_path.empty()) {
        fin.open(input_path);
        if (!fin) { std::cerr << "Cannot open: " << input_path << "\n"; return 1; }
        in = &fin;
    }
    std::ofstream fout;
    std::ostream* out = &std::cout;
    if (!output_path.empty()) {
        fout.open(output_path);
        if (!fout) { std::cerr << "Cannot open: " << output_path << "\n"; return 1; }
        out = &fout;
    }

    SensorManager                     sensor_mgr;
    InputParser                       parser(sensor_mgr);
    BlockingQueue<FullRawMeasurement> input_q(500'000);
    TOMHT                             tracker(cfg, sensor_mgr);

    OutputSerializer::write_header(*out);
    std::atomic<double> last_ts{0.0};

    // Thread 1: I/O reader
    std::thread io_thread([&]() {
        parser.stream_into(*in, input_q, verbose);
    });

    // Thread 2: Tracker
    std::thread tracker_thread([&]() {
        while (g_running.load(std::memory_order_relaxed)) {
            auto item = input_q.pop(50ms);
            if (!item) { if (input_q.is_stopped()) break; continue; }
            const SensorInfo& si = sensor_mgr.get_or_register(
                item->sensor_lat_deg, item->sensor_lon_deg, item->sensor_alt_m);
            Measurement meas = convert_measurement(*item, cfg.meas_noise_deg2);
            tracker.ingest(meas, si);
            last_ts.store(item->timestamp, std::memory_order_relaxed);
        }
    });

    // Thread 3: Output
    std::thread output_thread([&]() {
        while (g_running.load(std::memory_order_relaxed) || !input_q.is_stopped()) {
            std::this_thread::sleep_for(100ms);
            double ts = last_ts.load();
            auto snaps = tracker.get_outputs();
            if (!snaps.empty()) {
                OutputSerializer::write_batch(*out, ts, snaps);
                if (verbose)
                    for (const auto& s : snaps)
                        std::cerr << OutputSerializer::summary_line(ts, s) << "\n";
            }
        }
        double ts = last_ts.load();
        for (const auto& s : tracker.get_outputs())
            OutputSerializer::write_snapshot(*out, ts, s);
        out->flush();
    });

    io_thread.join();
    tracker_thread.join();
    g_running.store(false);
    output_thread.join();

    if (verbose)
        std::cerr << "[MTT] Done. Tracks=" << tracker.active_track_count() << "\n";
    return 0;
}
