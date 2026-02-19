#pragma once
#include <Eigen/Dense>
#include <cstdint>
#include <string>
#include <limits>
#include <cmath>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// WGS-84 Ellipsoid Constants
// ─────────────────────────────────────────────────────────────────────────────
namespace wgs84 {
    constexpr double a      = 6'378'137.0;              // semi-major axis [m]
    constexpr double f      = 1.0 / 298.257'223'563;    // flattening
    constexpr double b      = a * (1.0 - f);            // semi-minor axis [m]
    constexpr double e2     = 2.0*f - f*f;              // first eccentricity squared
    constexpr double ep2    = e2 / (1.0 - e2);         // second eccentricity squared
    constexpr double e      = std::sqrt(e2);
}  // namespace wgs84

// ─────────────────────────────────────────────────────────────────────────────
// Unit conversion
// ─────────────────────────────────────────────────────────────────────────────
constexpr double DEG2RAD = M_PI / 180.0;
constexpr double RAD2DEG = 180.0 / M_PI;

inline double deg2rad(double d) noexcept { return d * DEG2RAD; }
inline double rad2deg(double r) noexcept { return r * RAD2DEG; }

// ─────────────────────────────────────────────────────────────────────────────
// Eigen type aliases
// ─────────────────────────────────────────────────────────────────────────────
using Vec2d  = Eigen::Vector2d;
using Vec3d  = Eigen::Vector3d;
using Vec6d  = Eigen::Matrix<double, 6, 1>;
using Mat2d  = Eigen::Matrix2d;
using Mat3d  = Eigen::Matrix3d;
using Mat6d  = Eigen::Matrix<double, 6, 6>;
using Mat2x6 = Eigen::Matrix<double, 2, 6>;
using Mat6x2 = Eigen::Matrix<double, 6, 2>;
using MatXd  = Eigen::MatrixXd;
using VecXd  = Eigen::VectorXd;

// ─────────────────────────────────────────────────────────────────────────────
// Track ID
// ─────────────────────────────────────────────────────────────────────────────
using TrackID  = uint64_t;
using SensorID = std::string;   // keyed by "lat,lon,alt" string
constexpr TrackID INVALID_TRACK_ID = std::numeric_limits<TrackID>::max();

// ─────────────────────────────────────────────────────────────────────────────
// Sensor descriptor (pre-computed, cached)
// ─────────────────────────────────────────────────────────────────────────────
struct SensorInfo {
    SensorID  id;
    double    lat_deg, lon_deg, alt_m;   // geodetic (degrees / metres)
    Vec3d     ecef;                       // ECEF position [m]
    Vec3d     east_ecef;                  // ENU East  unit-vector in ECEF
    Vec3d     north_ecef;                 // ENU North unit-vector in ECEF
    Vec3d     up_ecef;                    // ENU Up    unit-vector in ECEF
};

// ─────────────────────────────────────────────────────────────────────────────
// Raw measurement (straight from the text file, in degrees)
// ─────────────────────────────────────────────────────────────────────────────
struct RawMeasurement {
    double    timestamp;      // seconds (or any monotonic time)
    double    az_deg;         // azimuth   [°]
    double    el_deg;         // elevation [°]
    SensorID  sensor_id;
};

// ─────────────────────────────────────────────────────────────────────────────
// Converted measurement (radians, ready for the filter)
// ─────────────────────────────────────────────────────────────────────────────
struct Measurement {
    double    timestamp;
    Vec2d     z;              // [azimuth_rad, elevation_rad]
    Mat2d     R;              // 2×2 measurement noise covariance [rad²]
    SensorID  sensor_id;
};

// ─────────────────────────────────────────────────────────────────────────────
// Track lifecycle states
// ─────────────────────────────────────────────────────────────────────────────
enum class TrackState : uint8_t {
    Tentative,    // newly born, not yet confirmed
    Confirmed,    // confirmed, actively output
    Coasting,     // no recent update – coasting on prediction only
    Deleted       // marked for return to object pool
};

// ─────────────────────────────────────────────────────────────────────────────
// MHT algorithm tuning knobs
// ─────────────────────────────────────────────────────────────────────────────
struct MHTConfig {
    // Gating
    double  gate_chi2_threshold   = 13.816; // chi2(2 DoF, 0.999)
    
    // Track management
    int     confirm_hits          = 3;      // consecutive updates to confirm
    int     deletion_misses       = 5;      // missed updates before deletion
    double  deletion_score_floor  = -15.0;  // log-likelihood score floor
    
    // TOMHT N-scan pruning
    int     n_scan_depth          = 5;
    int     k_best_hypotheses     = 10;
    
    // UKF sigma-point parameters
    double  ukf_alpha             = 1e-3;
    double  ukf_beta              = 2.0;
    double  ukf_kappa             = 0.0;
    
    // CWNA process noise power spectral density [m²/s³]
    double  psd_aircraft          = 1.0;    // smooth flight
    double  psd_bird              = 5.0;    // maneuvering biological
    // (clutter/pollen rejected by score drop, not separate PSD)
    double  psd_default           = 3.0;    // used until class. known
    
    // Measurement noise in degrees (from spec)
    double  meas_noise_deg2       = 0.01;   // σ² per axis [deg²]
    
    // Minimum sensors for triangulation
    int     min_sensors_for_init  = 2;
    
    // Initial velocity covariance [m²/s²]
    double  init_vel_cov          = 400.0;  // (20 m/s)²
};

}  // namespace mtt
