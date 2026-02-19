# MTT: Real-Time Multi-Target Tracker Design (C++)

## 1. Problem framing
MTT is a real-time **bearing-only** multi-target tracking system for 2-4 static cameras with 360° coverage. Input measurements are:

- `timestamp`
- `azimuth` (deg)
- `elevation` (deg)
- `sensor latitude` (deg, WGS-84)
- `sensor longitude` (deg, WGS-84)
- `sensor altitude` (m)

The environment can include hundreds of true airborne targets plus heavy clutter (dust/insects/pollen), creating a high false-alarm rate.

## 2. State and output definition
All tracks use a fully resolved 6D Cartesian state in ECEF:

`x = [px, py, pz, vx, vy, vz]`

where:
- position is meters in ECEF,
- velocity is meters/second in ECEF.

The CLI outputs **track ID** for every track together with ECEF position and ECEF velocity.

## 3. Multi-camera triangulation strategy
Because measurements are bearing-only, MTT resolves position by triangulating az/el rays from multiple cameras.

Pipeline per timestamp batch:
1. Convert each measurement to a sensor ray in ECEF (origin + direction).
2. Form cross-sensor ray pairs.
3. Compute closest approach for each pair of rays and midpoint candidate.
4. Reject weak intersections using:
   - max distance between paired rays,
   - az/el residual checks against both source measurements.
5. Cluster nearby candidates to produce one 3D point per target hypothesis.

This uses two or more camera measurements to produce resolved 3D positions robustly in clutter.

## 4. Track lifecycle and velocity estimation
The tracker performs:
- constant-velocity prediction in ECEF,
- nearest-neighbor association from predicted tracks to triangulated 3D candidates,
- alpha-beta style position/velocity update.

Velocity is estimated from finite-difference motion between consecutive associated 3D points and smoothed to suppress noise.

## 5. Clutter handling
Dense clutter is managed by a layered set of inexpensive filters:
- geometric cross-sensor consistency via triangulation pair separation,
- angular residual consistency checks,
- spatial clustering of candidate points,
- association gates and existence-probability decay for unobserved tracks,
- miss-streak based deletion.

## 6. Current repo mapping
- `include/mtt/types.hpp`: shared measurement/track/config types (ECEF 6D track state).
- `include/mtt/geo.hpp` + `src/geo.cpp`: WGS-84 + LOS conversions.
- `include/mtt/pmbm_tracker.hpp` + `src/pmbm_tracker.cpp`: tracker interface and triangulation-based implementation.
- `src/main.cpp`: CSV parsing, timestamp batching, tracker execution, and per-track ID output.

## 7. Next upgrades toward full PMBM
1. Replace nearest-neighbor association with PMBM global hypotheses and Murty k-best assignment.
2. Add UKF/cubature updates with explicit az/el measurement covariance usage.
3. Add triangulation confidence weighting and sensor-specific calibration/error models.
4. Add parallelized gating/association for high cardinality scenes.
5. Add OSPA/GOSPA evaluation and stress-test scenarios with dense clutter.
