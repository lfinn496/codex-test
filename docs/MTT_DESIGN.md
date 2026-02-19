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

## 2. Recommended state-of-the-art Bayesian core
For this regime, use a **PMBM filter** (Poisson Multi-Bernoulli Mixture):

- **Poisson intensity** models undetected targets and clutter-like births.
- **Multi-Bernoulli Mixture** models detected tracks with explicit existence probabilities.
- Naturally handles high clutter and unknown cardinality.
- Scales better than exhaustive MHT with appropriate pruning and capped hypotheses.

A practical implementation strategy:

1. Constant-velocity state in ECEF: `x = [px, py, pz, vx, vy, vz]`.
2. Nonlinear measurement model in az/el per sensor.
3. UKF or cubature updates per Bernoulli component.
4. Ellipsoidal gating in measurement space.
5. k-best assignment (Murty/Lagrangian) over gated candidates.
6. N-scan pruning + mixture reduction + track recycling to Poisson.

## 3. Coordinate and measurement model
### 3.1 Geodesy
- Convert each sensor pose (WGS-84 LLA) to ECEF.
- Represent LOS vectors by converting az/el in local ENU to ECEF direction vectors.

### 3.2 Measurement function
For track state `x`, sensor position `s`:

- Relative vector: `r = p - s`.
- Convert `r` to local ENU of sensor.
- Predicted azimuth/elevation:
  - `az = atan2(east, north)`
  - `el = atan2(up, sqrt(east^2 + north^2))`

The provided measurement covariance `R` is 2x2 in deg² and shared across all sensors.

## 4. Data association and clutter robustness
Use a layered strategy:

1. **Validation gating** with chi-square threshold in (az, el).
2. **Global assignment** with capped k-best hypotheses.
3. **Existence-probability logic** to suppress clutter-driven tracks.
4. **Track confirmation/deletion** rules (e.g., M/N confirmation, miss streak limits).
5. **Birth control**:
   - triangulation-consistent multi-sensor births favored,
   - isolated single-sensor births allowed with low prior existence.

## 5. Real-time engineering guidance
- Batch measurements by timestamp (or micro-batches in short windows).
- Parallelize per-track prediction and measurement likelihood computations.
- Bound compute with:
  - max hypotheses,
  - max tracks,
  - hard gating limits,
  - periodic pruning/merging.

For several hundred targets, this is mandatory for deterministic latency.

## 6. Repo implementation mapping
This repository includes a compile-ready scaffold:

- `include/mtt/types.hpp`: shared measurement/track/config types.
- `include/mtt/geo.hpp` + `src/geo.cpp`: WGS-84 and LOS conversions.
- `include/mtt/pmbm_tracker.hpp` + `src/pmbm_tracker.cpp`: PMBM-style tracker interface with extendable internal structure.
- `src/main.cpp`: CSV ingestion and tracker-driving CLI.

The current tracker implementation is intentionally lightweight and designed to be replaced incrementally with full PMBM internals while preserving interfaces.

## 7. Near-term roadmap
1. Replace placeholder update with full Bernoulli + global-hypothesis update.
2. Add nonlinear measurement Jacobian/UKF and angular wrap-safe innovations.
3. Add multi-sensor triangulation-based birth initialization.
4. Add configurable clutter map/intensity over az-el domain.
5. Add evaluation scripts for OSPA/GOSPA and CLEAR MOT-like metrics on logged truth.
