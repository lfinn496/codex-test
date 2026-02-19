// mtt/tests/test_main.cpp
// Lightweight unit tests. Run: ./mtt_test  (exits 0 on pass, nonzero on fail)

#include "mtt/types.hpp"
#include "mtt/geodesy.hpp"
#include "mtt/observation_model.hpp"
#include "mtt/ukf.hpp"
#include "mtt/ws3d.hpp"
#include "mtt/assignment.hpp"
#include "mtt/gating.hpp"
#include "mtt/track.hpp"
#include "mtt/input_parser.hpp"
#include "mtt/sensor_manager.hpp"

#include <cmath>
#include <iostream>
#include <sstream>
#include <cassert>
#include <stdexcept>

using namespace mtt;

// ─── Test harness ─────────────────────────────────────────────────────────────
static int pass_count = 0;
static int fail_count = 0;

#define ASSERT_NEAR(a, b, tol, msg) do {                                \
    double _a = (a), _b = (b);                                          \
    if (std::abs(_a - _b) > (tol)) {                                    \
        std::cerr << "[FAIL] " << (msg)                                 \
                  << "  got=" << _a << " expected=" << _b               \
                  << " diff=" << std::abs(_a - _b) << "\n";             \
        ++fail_count;                                                   \
    } else {                                                            \
        std::cout << "[PASS] " << (msg) << "\n";                        \
        ++pass_count;                                                   \
    }                                                                   \
} while (0)

#define ASSERT_TRUE(cond, msg) do {                                     \
    if (!(cond)) {                                                      \
        std::cerr << "[FAIL] " << (msg) << "\n";                        \
        ++fail_count;                                                   \
    } else {                                                            \
        std::cout << "[PASS] " << (msg) << "\n";                        \
        ++pass_count;                                                   \
    }                                                                   \
} while (0)

// ─── Helper: build a BearingRay from a SensorInfo + true target ECEF ─────────
static BearingRay make_ray(const SensorInfo& s, const Vec3d& target_ecef,
                           double sigma2 = (0.1 * DEG2RAD) * (0.1 * DEG2RAD))
{
    Vec6d st = Vec6d::Zero();
    st.head<3>() = target_ecef;
    Vec2d z = bearing_observation(st, s);
    BearingRay r;
    r.sensor_ecef = s.ecef;
    r.az_rad      = z(0);
    r.el_rad      = z(1);
    r.sigma2_az   = sigma2;
    r.sigma2_el   = sigma2;
    r.east        = s.east_ecef;
    r.north       = s.north_ecef;
    r.up          = s.up_ecef;
    return r;
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. Geodesy: round-trip LLA → ECEF → LLA
// ─────────────────────────────────────────────────────────────────────────────
void test_geodesy_roundtrip() {
    const double lat0 = 51.4779, lon0 = -0.0015, alt0 = 46.0;
    Vec3d ecef = geodesy::lla_to_ecef(lat0, lon0, alt0);
    double lat1, lon1, alt1;
    geodesy::ecef_to_lla(ecef, lat1, lon1, alt1);
    ASSERT_NEAR(lat1, lat0, 1e-6, "geodesy roundtrip lat (deg)");
    ASSERT_NEAR(lon1, lon0, 1e-6, "geodesy roundtrip lon (deg)");
    ASSERT_NEAR(alt1, alt0, 0.01, "geodesy roundtrip alt (m)");
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. ENU axes must be orthonormal
// ─────────────────────────────────────────────────────────────────────────────
void test_enu_axes_orthonormal() {
    Vec3d e, n, u;
    geodesy::compute_enu_axes(48.8566, 2.3522, e, n, u);
    ASSERT_NEAR(e.norm(), 1.0, 1e-12, "ENU East unit norm");
    ASSERT_NEAR(n.norm(), 1.0, 1e-12, "ENU North unit norm");
    ASSERT_NEAR(u.norm(), 1.0, 1e-12, "ENU Up unit norm");
    ASSERT_NEAR(e.dot(n), 0.0, 1e-12, "ENU E.N orthogonal");
    ASSERT_NEAR(e.dot(u), 0.0, 1e-12, "ENU E.U orthogonal");
    ASSERT_NEAR(n.dot(u), 0.0, 1e-12, "ENU N.U orthogonal");
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Observation model: known position → known bearing
// ─────────────────────────────────────────────────────────────────────────────
void test_observation_model() {
    const SensorInfo s = geodesy::make_sensor_info("test", 0.0, 0.0, 0.0);
    Vec6d state = Vec6d::Zero();
    state.head<3>() = s.ecef + s.north_ecef * 1000.0 + s.up_ecef * 100.0;
    Vec2d z = bearing_observation(state, s);
    ASSERT_NEAR(z(0), 0.0, 1e-3, "obs model: due north azimuth ~= 0");
    const double expected_el = std::atan2(100.0, 1000.0);
    ASSERT_NEAR(z(1), expected_el, 1e-4, "obs model: elevation atan2(100,1000)");
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. Measurement conversion: degrees^2 -> radians^2
// ─────────────────────────────────────────────────────────────────────────────
void test_measurement_conversion() {
    RawMeasurement raw;
    raw.timestamp = 1.0;
    raw.az_deg    = 45.0;
    raw.el_deg    = 30.0;
    raw.sensor_id = "s0";
    Measurement m = convert_measurement(raw, 1.0);
    ASSERT_NEAR(m.z(0), deg2rad(45.0), 1e-12, "meas conv: azimuth rad");
    ASSERT_NEAR(m.z(1), deg2rad(30.0), 1e-12, "meas conv: elevation rad");
    const double expected_r = DEG2RAD * DEG2RAD;
    ASSERT_NEAR(m.R(0,0), expected_r, 1e-20, "meas conv: R diagonal deg->rad");
    ASSERT_NEAR(m.R(1,1), expected_r, 1e-20, "meas conv: R off-diagonal zero");
    ASSERT_NEAR(m.R(0,1), 0.0,        1e-20, "meas conv: R cross zero");
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. CWNA matrices: symmetry, positive-definiteness
// ─────────────────────────────────────────────────────────────────────────────
void test_cwna_matrices() {
    const double dt = 0.033;
    const double q  = 2.0;
    Mat6d F = UKF::build_F(dt);
    Mat6d Q = UKF::build_Q(dt, q);
    ASSERT_NEAR(F(0,3), dt,  1e-12, "CWNA F(0,3)=dt");
    ASSERT_NEAR(F(3,0), 0.0, 1e-12, "CWNA F(3,0)=0");
    const Mat6d diff = Q - Q.transpose();
    ASSERT_NEAR(diff.norm(), 0.0, 1e-12, "CWNA Q symmetric");
    Eigen::LLT<Mat6d> llt(Q);
    ASSERT_TRUE(llt.info() == Eigen::Success, "CWNA Q positive-definite");
}

// ─────────────────────────────────────────────────────────────────────────────
// 6. UKF predict+update stays finite and P stays PD
// ─────────────────────────────────────────────────────────────────────────────
void test_ukf_predict_update() {
    MHTConfig cfg;
    UKF filter(cfg);
    UKFState s;
    SensorInfo sensor = geodesy::make_sensor_info("s", 51.4779, -0.0015, 46.0);
    s.x.head<3>() = sensor.ecef + sensor.north_ecef * 10000.0 + sensor.up_ecef * 1000.0;
    s.P = Mat6d::Identity() * 1000.0;
    s.P.bottomRightCorner<3,3>() *= 100.0;
    s.timestamp = 0.0;
    filter.predict(s, 0.1);
    Vec2d z_true = bearing_observation(s.x, sensor);
    Measurement meas;
    meas.timestamp = 0.1;
    meas.sensor_id = sensor.id;
    meas.z = z_true;
    meas.z(0) += 1e-4;
    meas.z(1) += 1e-4;
    meas.R = Mat2d::Identity() * (0.01 * DEG2RAD * DEG2RAD);
    double maha2 = filter.update(s, meas, sensor);
    ASSERT_TRUE(std::isfinite(maha2),     "UKF update: maha2 finite");
    ASSERT_TRUE(std::isfinite(s.x.sum()), "UKF update: state finite");
    Eigen::LLT<Mat6d> llt(s.P);
    ASSERT_TRUE(llt.info() == Eigen::Success, "UKF update: P PD after update");
}

// ─────────────────────────────────────────────────────────────────────────────
// 7. Hungarian algorithm: 3x3 cost matrix
// ─────────────────────────────────────────────────────────────────────────────
void test_hungarian_3x3() {
    HungarianSolver::CostMatrix cost = {
        {9, 2, 7},
        {3, 6, 2},
        {6, 3, 4}
    };
    auto asgn = HungarianSolver::solve(cost);
    ASSERT_TRUE(asgn.size() == 3, "Hungarian: result size 3");
    std::vector<bool> used(3, false);
    bool unique = true;
    for (int a : asgn) {
        if (a < 0 || a >= 3 || used[a]) { unique = false; break; }
        used[a] = true;
    }
    ASSERT_TRUE(unique, "Hungarian: unique assignment");
}

// ─────────────────────────────────────────────────────────────────────────────
// 8. WS3D 2-sensor: well-conditioned geometry (sensors not co-linear w/ target)
//
// Baseline: s1 at (0,0), s2 at (0.009 deg lat, 0) => ~1 km North.
// Target:   1 km East of s1 midpoint, 500 m above ground.
// Azimuths: ~90 deg from s1 (due East), ~135 deg from s2 (SE), non-degenerate.
// ─────────────────────────────────────────────────────────────────────────────
void test_ws3d_2sensor() {
    SensorInfo s1 = geodesy::make_sensor_info("s1", 0.0,   0.0,   0.0);
    SensorInfo s2 = geodesy::make_sensor_info("s2", 0.009, 0.0,   0.0); // ~1 km North

    // Target: 1 km East of s1, 500 m altitude
    const Vec3d target = s1.ecef
                       + s1.east_ecef  * 1000.0
                       + s1.up_ecef    *  500.0;

    WS3DResult res = WS3DTriangulator::triangulate({make_ray(s1, target),
                                                    make_ray(s2, target)});
    ASSERT_TRUE(res.valid, "WS3D 2-sensor: valid result");

    const double err = (res.position_ecef - target).norm();
    if (err >= 50.0)
        std::cerr << "  WS3D 2-sensor error=" << err << " m\n";
    ASSERT_TRUE(err < 50.0, "WS3D 2-sensor: position error < 50 m");
}

// ─────────────────────────────────────────────────────────────────────────────
// 9. WS3D 4-sensor: sensors at corners of a 2 km square
// ─────────────────────────────────────────────────────────────────────────────
void test_ws3d_4sensor() {
    SensorInfo s1 = geodesy::make_sensor_info("s1", 0.0,   0.0,   0.0);
    SensorInfo s2 = geodesy::make_sensor_info("s2", 0.0,   0.018, 0.0);
    SensorInfo s3 = geodesy::make_sensor_info("s3", 0.018, 0.0,   0.0);
    SensorInfo s4 = geodesy::make_sensor_info("s4", 0.018, 0.018, 0.0);

    const Vec3d centre = (s1.ecef + s2.ecef + s3.ecef + s4.ecef) * 0.25
                       + s1.up_ecef * 1000.0;

    WS3DResult res = WS3DTriangulator::triangulate({make_ray(s1, centre),
                                                    make_ray(s2, centre),
                                                    make_ray(s3, centre),
                                                    make_ray(s4, centre)});
    ASSERT_TRUE(res.valid, "WS3D 4-sensor: valid result");
    const double err = (res.position_ecef - centre).norm();
    if (err >= 100.0)
        std::cerr << "  WS3D 4-sensor error=" << err << " m\n";
    ASSERT_TRUE(err < 100.0, "WS3D 4-sensor: position error < 100 m");
}

// ─────────────────────────────────────────────────────────────────────────────
// 10. Input parser: well-formed line
// ─────────────────────────────────────────────────────────────────────────────
void test_input_parser() {
    SensorManager sm;
    InputParser parser(sm);
    FullRawMeasurement m;
    bool ok = parser.parse_line("1000.5 45.123 12.456 51.477 -0.001 100.0", m);
    ASSERT_TRUE(ok,                          "InputParser: valid line parsed");
    ASSERT_NEAR(m.timestamp, 1000.5, 1e-9,  "InputParser: timestamp");
    ASSERT_NEAR(m.az_deg,    45.123, 1e-9,  "InputParser: azimuth");
    ASSERT_NEAR(m.el_deg,    12.456, 1e-9,  "InputParser: elevation");
}

// ─────────────────────────────────────────────────────────────────────────────
// 11. Coarse gating: within-gate and out-of-gate
// ─────────────────────────────────────────────────────────────────────────────
void test_gating() {
    MHTConfig cfg;
    UKF filter(cfg);
    CoarseGate gate(cfg.gate_chi2_threshold);

    SensorInfo sensor = geodesy::make_sensor_info("s", 51.4779, -0.0015, 46.0);
    Track t;
    t.ukf.x.head<3>() = sensor.ecef + sensor.north_ecef * 5000.0 + sensor.up_ecef * 500.0;
    t.ukf.P = Mat6d::Identity() * 1e6;
    t.ukf.timestamp = 0.0;

    Vec2d z_pred = bearing_observation(t.ukf.x, sensor);

    Measurement m_in;
    m_in.z = z_pred;
    m_in.R = Mat2d::Identity() * (0.01 * DEG2RAD * DEG2RAD);
    m_in.timestamp = 0.0;
    m_in.sensor_id = sensor.id;

    std::vector<Track*> tracks = {&t};
    auto gated_in = gate.filter(m_in, sensor, tracks, filter);
    ASSERT_TRUE(!gated_in.empty(), "Gating: predicted meas passes gate");

    Measurement m_out = m_in;
    m_out.z(0) += M_PI / 2.0;
    auto gated_out = gate.filter(m_out, sensor, tracks, filter);
    ASSERT_TRUE(gated_out.empty(), "Gating: distant meas fails gate");
}

// ─────────────────────────────────────────────────────────────────────────────
int main() {
    std::cout << "=== MTT Unit Tests ===\n\n";
    test_geodesy_roundtrip();
    test_enu_axes_orthonormal();
    test_observation_model();
    test_measurement_conversion();
    test_cwna_matrices();
    test_ukf_predict_update();
    test_hungarian_3x3();
    test_ws3d_2sensor();
    test_ws3d_4sensor();
    test_input_parser();
    test_gating();
    std::cout << "\n=== Results: "
              << pass_count << " passed, "
              << fail_count << " failed ===\n";
    return fail_count > 0 ? 1 : 0;
}
