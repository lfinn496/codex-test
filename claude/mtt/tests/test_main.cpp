// mtt/tests/test_main.cpp
// Lightweight unit tests – no external framework needed.
// Run with: ./mtt_test  (exits 0 on pass, nonzero on failure)

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

#define ASSERT_NEAR(a, b, tol, msg) do {                         \
    double _a = (a), _b = (b);                                   \
    if (std::abs(_a - _b) > (tol)) {                             \
        std::cerr << "[FAIL] " << (msg)                          \
                  << "  got=" << _a << " expected=" << _b        \
                  << " diff=" << std::abs(_a-_b) << "\n";        \
        ++fail_count;                                             \
    } else {                                                      \
        std::cout << "[PASS] " << (msg) << "\n";                 \
        ++pass_count;                                             \
    }                                                             \
} while(0)

#define ASSERT_TRUE(cond, msg) do {                              \
    if (!(cond)) {                                               \
        std::cerr << "[FAIL] " << (msg) << "\n";                 \
        ++fail_count;                                            \
    } else {                                                     \
        std::cout << "[PASS] " << (msg) << "\n";                 \
        ++pass_count;                                            \
    }                                                            \
} while(0)

// ─── 1. Geodesy: round-trip LLA → ECEF → LLA ─────────────────────────────────
void test_geodesy_roundtrip() {
    const double lat0 =  51.4779,  lon0 = -0.0015, alt0 = 46.0;  // Greenwich
    Vec3d ecef = geodesy::lla_to_ecef(lat0, lon0, alt0);

    double lat1, lon1, alt1;
    geodesy::ecef_to_lla(ecef, lat1, lon1, alt1);

    ASSERT_NEAR(lat1, lat0, 1e-6, "geodesy roundtrip lat (deg)");
    ASSERT_NEAR(lon1, lon0, 1e-6, "geodesy roundtrip lon (deg)");
    ASSERT_NEAR(alt1, alt0, 0.01, "geodesy roundtrip alt (m)");
}

// ENU axes must be orthonormal
void test_enu_axes_orthonormal() {
    Vec3d e, n, u;
    geodesy::compute_enu_axes(48.8566, 2.3522, e, n, u);  // Paris
    ASSERT_NEAR(e.norm(), 1.0, 1e-12, "ENU East unit norm");
    ASSERT_NEAR(n.norm(), 1.0, 1e-12, "ENU North unit norm");
    ASSERT_NEAR(u.norm(), 1.0, 1e-12, "ENU Up unit norm");
    ASSERT_NEAR(e.dot(n), 0.0, 1e-12, "ENU E·N orthogonal");
    ASSERT_NEAR(e.dot(u), 0.0, 1e-12, "ENU E·U orthogonal");
    ASSERT_NEAR(n.dot(u), 0.0, 1e-12, "ENU N·U orthogonal");
}

// ─── 2. Observation model: known position should yield known bearing ───────────
void test_observation_model() {
    // Sensor at origin of a local ENU frame, target due north at 1000 m, 100 m up
    const double lat = 0.0, lon = 0.0, alt = 0.0;
    SensorInfo s = geodesy::make_sensor_info("test", lat, lon, alt);

    // Target: sensor ECEF + 1000*north + 100*up
    Vec6d state = Vec6d::Zero();
    state.head<3>() = s.ecef + s.north_ecef * 1000.0 + s.up_ecef * 100.0;

    Vec2d z = bearing_observation(state, s);
    // Azimuth: target is due north → az ≈ 0
    ASSERT_NEAR(z(0), 0.0, 1e-3, "obs model: due north azimuth ≈ 0");
    // Elevation: tan(el) = 100/1000 → el ≈ 5.71°
    const double expected_el = std::atan2(100.0, 1000.0);
    ASSERT_NEAR(z(1), expected_el, 1e-4, "obs model: elevation atan2(100,1000)");
}

// ─── 3. Measurement conversion: degrees² → radians² ──────────────────────────
void test_measurement_conversion() {
    RawMeasurement raw;
    raw.timestamp = 1.0;
    raw.az_deg    = 45.0;
    raw.el_deg    = 30.0;
    raw.sensor_id = "s0";

    Measurement m = convert_measurement(raw, 1.0);   // σ²=1 deg²
    ASSERT_NEAR(m.z(0), deg2rad(45.0), 1e-12, "meas conv: azimuth rad");
    ASSERT_NEAR(m.z(1), deg2rad(30.0), 1e-12, "meas conv: elevation rad");
    const double expected_r = DEG2RAD * DEG2RAD;
    ASSERT_NEAR(m.R(0,0), expected_r, 1e-20, "meas conv: R diagonal deg→rad");
    ASSERT_NEAR(m.R(1,1), expected_r, 1e-20, "meas conv: R off-diagonal zero");
    ASSERT_NEAR(m.R(0,1), 0.0,        1e-20, "meas conv: R cross zero");
}

// ─── 4. CWNA matrices: symmetry, positive-definiteness ───────────────────────
void test_cwna_matrices() {
    const double dt = 0.033;   // ~30 fps
    const double q  = 2.0;
    Mat6d F = UKF::build_F(dt);
    Mat6d Q = UKF::build_Q(dt, q);

    // F top-right should be dt*I
    ASSERT_NEAR(F(0,3), dt, 1e-12, "CWNA F(0,3)=dt");
    ASSERT_NEAR(F(3,0), 0.0, 1e-12, "CWNA F(3,0)=0");

    // Q must be symmetric PD
    const Mat6d diff = Q - Q.transpose();
    ASSERT_NEAR(diff.norm(), 0.0, 1e-12, "CWNA Q symmetric");

    Eigen::LLT<Mat6d> llt(Q);
    ASSERT_TRUE(llt.info() == Eigen::Success, "CWNA Q positive-definite");
}

// ─── 5. UKF predict+update doesn't blow up ────────────────────────────────────
void test_ukf_predict_update() {
    MHTConfig cfg;
    UKF filter(cfg);

    UKFState s;
    // Place target 10 km north of Greenwich at 1000 m altitude
    SensorInfo sensor = geodesy::make_sensor_info("s", 51.4779, -0.0015, 46.0);
    s.x.head<3>() = sensor.ecef + sensor.north_ecef * 10'000.0 + sensor.up_ecef * 1000.0;
    s.P = Mat6d::Identity() * 1000.0;
    s.P.bottomRightCorner<3,3>() *= 100.0;
    s.timestamp = 0.0;

    // Predict 0.1s forward
    filter.predict(s, 0.1);

    // Generate a synthetic measurement
    Vec2d z_true = bearing_observation(s.x, sensor);

    // Build a measurement with small noise
    Measurement meas;
    meas.timestamp = 0.1;
    meas.sensor_id = sensor.id;
    meas.z = z_true;
meas.z(0) += 1e-4; meas.z(1) += 1e-4;  // tiny perturbation
    meas.R = Mat2d::Identity() * (0.01 * DEG2RAD * DEG2RAD);

    double maha2 = filter.update(s, meas, sensor);

    ASSERT_TRUE(std::isfinite(maha2),    "UKF update: maha2 finite");
    ASSERT_TRUE(std::isfinite(s.x.sum()),"UKF update: state finite");

    // Covariance should still be PD
    Eigen::LLT<Mat6d> llt(s.P);
    ASSERT_TRUE(llt.info() == Eigen::Success, "UKF update: P PD after update");
}

// ─── 6. Hungarian algorithm: simple 3×3 cost ─────────────────────────────────
void test_hungarian_3x3() {
    // Classic example: optimal cost = 11 (0→2=6, 1→0=2, 2→1=3)
    HungarianSolver::CostMatrix cost = {
        {9, 2, 7},
        {3, 6, 2},
        {6, 3, 4}
    };
    auto asgn = HungarianSolver::solve(cost);
    // Each row assigned to a unique column
    ASSERT_TRUE(asgn.size() == 3, "Hungarian: result size 3");
    // Verify uniqueness
    std::vector<bool> used(3, false);
    bool unique = true;
    for (int a : asgn) {
        if (a < 0 || a >= 3 || used[a]) { unique = false; break; }
        used[a] = true;
    }
    ASSERT_TRUE(unique, "Hungarian: unique assignment");
}

// ─── 7. WS3D triangulation: two sensors, known geometry ──────────────────────
void test_ws3d_triangulation() {
    // Two sensors 1 km apart (east-west), target at their midpoint + 500 m up
    SensorInfo s1 = geodesy::make_sensor_info("s1",  0.0,  0.0,  0.0);
    SensorInfo s2 = geodesy::make_sensor_info("s2",  0.0,  0.009, 0.0);  // ~1 km east

    // True target in ECEF
    const Vec3d true_pos = (s1.ecef + s2.ecef) * 0.5 + s1.up_ecef * 500.0;

    // Compute bearings from each sensor
    Vec6d state = Vec6d::Zero();
    state.head<3>() = true_pos;
    Vec2d z1 = bearing_observation(state, s1);
    Vec2d z2 = bearing_observation(state, s2);

    const double sigma2 = (0.1 * DEG2RAD) * (0.1 * DEG2RAD);
    BearingRay r1{ s1.ecef, z1(0), z1(1), sigma2, sigma2, s1.east_ecef, s1.north_ecef, s1.up_ecef };
    BearingRay r2{ s2.ecef, z2(0), z2(1), sigma2, sigma2, s2.east_ecef, s2.north_ecef, s2.up_ecef };

    WS3DResult res = WS3DTriangulator::triangulate({r1, r2});
    ASSERT_TRUE(res.valid, "WS3D: valid result");

    const double err = (res.position_ecef - true_pos).norm();
    // Expect within ~200 m for 0.1° noise at ~500 m range
    ASSERT_TRUE(err < 500.0, "WS3D: position error < 500 m");
}

// ─── 8. Input parser: well-formed line ───────────────────────────────────────
void test_input_parser() {
    SensorManager sm;
    InputParser parser(sm);
    FullRawMeasurement m;
    bool ok = parser.parse_line("1000.5 45.123 12.456 51.477 -0.001 100.0", m);
    ASSERT_TRUE(ok, "InputParser: valid line parsed");
    ASSERT_NEAR(m.timestamp, 1000.5, 1e-9, "InputParser: timestamp");
    ASSERT_NEAR(m.az_deg,    45.123, 1e-9, "InputParser: azimuth");
    ASSERT_NEAR(m.el_deg,    12.456, 1e-9, "InputParser: elevation");
}

// ─── 9. Coarse gating: within-gate and out-of-gate ────────────────────────────
void test_gating() {
    MHTConfig cfg;
    UKF filter(cfg);
    CoarseGate gate(cfg.gate_chi2_threshold);

    // Build a track at a known position
    SensorInfo sensor = geodesy::make_sensor_info("s", 51.4779, -0.0015, 46.0);
    Track t;
    t.ukf.x.head<3>() = sensor.ecef + sensor.north_ecef * 5000.0 + sensor.up_ecef * 500.0;
    t.ukf.P = Mat6d::Identity() * 1e6;
    t.ukf.timestamp = 0.0;

    // Measurement at predicted location (should pass gate)
    Vec2d z_pred = bearing_observation(t.ukf.x, sensor);
    Measurement m_in;
    m_in.z = z_pred;
    m_in.R = Mat2d::Identity() * (0.01 * DEG2RAD * DEG2RAD);
    m_in.timestamp = 0.0;
    m_in.sensor_id = sensor.id;

    std::vector<Track*> tracks = {&t};
    auto gated_in = gate.filter(m_in, sensor, tracks, filter);
    ASSERT_TRUE(!gated_in.empty(), "Gating: predicted meas passes gate");

    // Measurement 90° away (should fail gate)
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
    test_ws3d_triangulation();
    test_input_parser();
    test_gating();

    std::cout << "\n=== Results: "
              << pass_count << " passed, "
              << fail_count << " failed ===\n";

    return fail_count > 0 ? 1 : 0;
}
