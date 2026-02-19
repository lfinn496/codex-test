#pragma once
#include "types.hpp"
#include <vector>
#include <cmath>
#include <Eigen/QR>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// Single bearing ray used by WS3D
// ─────────────────────────────────────────────────────────────────────────────
struct BearingRay {
    Vec3d  sensor_ecef;   // sensor position [m]
    double az_rad;        // azimuth   [rad]
    double el_rad;        // elevation [rad]
    double sigma2_az;     // azimuth   variance [rad²]
    double sigma2_el;     // elevation variance [rad²]
    Vec3d  east;          // sensor ENU east  in ECEF
    Vec3d  north;         // sensor ENU north in ECEF
    Vec3d  up;            // sensor ENU up    in ECEF
};

// ─────────────────────────────────────────────────────────────────────────────
// Result of WS3D triangulation
// ─────────────────────────────────────────────────────────────────────────────
struct WS3DResult {
    Vec3d   position_ecef;  // triangulated 3D position [m]
    Mat3d   cov_pos;        // 3×3 position covariance [m²]
    double  quality;        // trace(cov_pos)
    bool    valid = false;
};

// ─────────────────────────────────────────────────────────────────────────────
// Weighted 3-D bearing-only triangulation.
//
// For each ray the 3-D unit direction vector in ENU is:
//   d = (sin(az)*cos(el),  cos(az)*cos(el),  sin(el))
//
// The target position x minimises the sum of squared perpendicular distances
// from each ray:
//   min  Σ_i  w_i * ||(I - d_i d_i^T)(x - q_i)||²
//
// Normal equations:
//   A x = b,   A = Σ w_i (I - d_i d_i^T),   b = Σ w_i (I - d_i d_i^T) q_i
//
// All vectors are expressed in the ENU frame of sensor[0] for numerical
// stability (offsets are O(km) rather than O(6000 km)), then rotated back
// to ECEF.
//
// This formulation is never rank-deficient as long as the rays are not all
// parallel, and subsumes the original two-step Stansfield approach without
// its azimuth-degeneracy failure mode.
// ─────────────────────────────────────────────────────────────────────────────
class WS3DTriangulator {
public:
    static constexpr int MIN_RAYS = 2;

    static WS3DResult triangulate(const std::vector<BearingRay>& rays) {
        const int m = static_cast<int>(rays.size());
        if (m < MIN_RAYS) return {};

        // Reference sensor defines the local ENU origin
        const BearingRay& ref = rays[0];

        // Build normal-equation matrices in ENU
        Mat3d A = Mat3d::Zero();
        Vec3d b = Vec3d::Zero();

        for (int i = 0; i < m; ++i) {
            // 3-D unit direction in ENU
            const double cos_el = std::cos(rays[i].el_rad);
            const double sin_el = std::sin(rays[i].el_rad);
            const double cos_az = std::cos(rays[i].az_rad);
            const double sin_az = std::sin(rays[i].az_rad);
            Vec3d d;
            d << sin_az * cos_el,
                 cos_az * cos_el,
                 sin_el;

            // Sensor position in reference ENU frame
            const Vec3d dr = rays[i].sensor_ecef - ref.sensor_ecef;
            Vec3d q;
            q << ref.east.dot(dr),
                 ref.north.dot(dr),
                 ref.up.dot(dr);

            // Weight: 1 / (σ²_az + σ²_el) – penalises noisier sensors
            const double w = 1.0 / (rays[i].sigma2_az + rays[i].sigma2_el + 1e-30);

            // I - d d^T
            Mat3d P = Mat3d::Identity();
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    P(r, c) -= d(r) * d(c);

            A = A + P * w;
            b = b + P * q * w;
        }

        // Solve A x = b in ENU; add ridge for numerical stability
        A = A + Mat3d::Identity() * 1e-6;
        const Vec3d x_enu = A.ldlt().solve(b);

        // Rotate ENU result back to ECEF
        const Vec3d pos_ecef = ref.sensor_ecef
                             + ref.east  * x_enu(0)
                             + ref.north * x_enu(1)
                             + ref.up    * x_enu(2);

        // Approximate position covariance: A^{-1} in ENU, then rotate to ECEF
        const Mat3d cov_enu = A.inverse();

        // Rotation matrix R: columns are east, north, up unit vectors
        Mat3d R;
        for (int r = 0; r < 3; ++r) {
            R(r, 0) = ref.east(r);
            R(r, 1) = ref.north(r);
            R(r, 2) = ref.up(r);
        }
        const Mat3d cov_ecef = R * cov_enu * R.transpose();

        WS3DResult res;
        res.position_ecef = pos_ecef;
        res.cov_pos       = cov_ecef;
        res.quality       = cov_ecef.trace();
        res.valid         = true;
        return res;
    }
};

}  // namespace mtt
