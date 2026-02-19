#pragma once
#include "types.hpp"
#include "track.hpp"
#include <ostream>
#include <iomanip>
#include <string>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// Output columns (space-separated):
//
//   timestamp  track_id  state  Px  Py  Pz  Vx  Vy  Vz
//
// where P = ECEF position [m], V = ECEF velocity [m/s].
// state: 0=Tentative 1=Confirmed 2=Coasting
// ─────────────────────────────────────────────────────────────────────────────
class OutputSerializer {
public:
    // Write CSV header
    static void write_header(std::ostream& os) {
        os << "# timestamp,track_id,state,"
           << "px_m,py_m,pz_m,"
           << "vx_mps,vy_mps,vz_mps\n";
    }

    static void write_snapshot(std::ostream& os,
                                double timestamp,
                                const Track::Snapshot& s) {
        os << std::fixed << std::setprecision(6)
           << timestamp                       << ','
           << s.id                             << ','
           << static_cast<int>(s.state)        << ','
           << std::setprecision(3)
           << s.state_vec(0)  << ','           // Px
           << s.state_vec(1)  << ','           // Py
           << s.state_vec(2)  << ','           // Pz
           << std::setprecision(4)
           << s.state_vec(3)  << ','           // Vx
           << s.state_vec(4)  << ','           // Vy
           << s.state_vec(5)  << '\n';         // Vz
    }

    static void write_batch(std::ostream& os,
                             double timestamp,
                             const std::vector<Track::Snapshot>& snaps) {
        for (const auto& s : snaps)
            write_snapshot(os, timestamp, s);
        os.flush();
    }

    // Human-readable summary line for console monitoring
    static std::string summary_line(double timestamp,
                                    const Track::Snapshot& s) {
        char buf[256];
        const Vec3d p = s.state_vec.head<3>();
        const Vec3d v = s.state_vec.tail<3>();

        // Convert ECEF position back to LLA for readability
        double lat, lon, alt;
        // (Inline Bowring – avoids including geodesy header here)
        const double x = p(0), y = p(1), z = p(2);
        lon = std::atan2(y, x) * RAD2DEG;
        const double p2 = std::hypot(x, y);
        double lat_r = std::atan2(z, p2 * (1.0 - wgs84::e2));
        for (int i = 0; i < 5; ++i) {
            double N = wgs84::a / std::sqrt(1.0 - wgs84::e2*std::sin(lat_r)*std::sin(lat_r));
            lat_r = std::atan2(z + wgs84::e2 * N * std::sin(lat_r), p2);
        }
        lat = lat_r * RAD2DEG;
        const double N = wgs84::a / std::sqrt(1.0 - wgs84::e2*std::sin(lat_r)*std::sin(lat_r));
        alt = (std::cos(lat_r) > 1e-9) ? p2/std::cos(lat_r) - N : 0.0;

        const double speed_mps = v.norm();
        snprintf(buf, sizeof(buf),
            "t=%.3f TID=%5llu state=%d lat=%.5f lon=%.5f alt=%.1fm spd=%.1fm/s",
            timestamp, (unsigned long long)s.id,
            (int)s.state, lat, lon, alt, speed_mps);
        return buf;
    }
};

}  // namespace mtt
