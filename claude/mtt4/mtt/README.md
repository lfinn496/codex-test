# MTT — Real-Time Multi-Target Tracker

C++17 implementation of the bearing-only multi-target tracking system specified
in *"Architectural and Mathematical Blueprint for the MTT Real-Time
Multi-Target Tracking System."*

---

## Architecture

```
Input text stream
       │
       ▼
 ┌─────────────┐     ┌─────────────────────────────────────────────────────────┐
 │ InputParser │────▶│ BlockingQueue<FullRawMeasurement>  (500k capacity)       │
 │ (I/O Thread)│     └──────────────────────┬──────────────────────────────────┘
 └─────────────┘                            │
                                            ▼
                              ┌─────────────────────────┐
                              │     Tracker Thread       │
                              │                          │
                              │  SensorManager           │
                              │  (LLA→ECEF cache)        │
                              │         │                │
                              │         ▼                │
                              │  convert_measurement()   │
                              │  (deg → rad, R matrix)   │
                              │         │                │
                              │         ▼                │
                              │  ┌─────────────────────┐ │
                              │  │     TOMHT Engine     │ │
                              │  │                      │ │
                              │  │  Ray buffer          │ │
                              │  │       ↓              │ │
                              │  │  WS3D triangulation  │ │
                              │  │  (track birth)       │ │
                              │  │       ↓              │ │
                              │  │  CoarseGate          │ │
                              │  │  (Mahalanobis χ²)    │ │
                              │  │       ↓              │ │
                              │  │  Hungarian assign.   │ │
                              │  │       ↓              │ │
                              │  │  UKF sequential upd. │ │
                              │  │  (CWNA predict +     │ │
                              │  │   sigma-pt update)   │ │
                              │  │       ↓              │ │
                              │  │  Track lifecycle mgr │ │
                              │  │  (confirm/delete/    │ │
                              │  │   pool return)       │ │
                              │  └─────────────────────┘ │
                              └─────────────────────────┘
                                            │
                                            ▼
                              ┌─────────────────────────┐
                              │     Output Thread        │
                              │  OutputSerializer        │
                              │  (track_id, 6D ECEF)     │
                              └─────────────────────────┘
                                            │
                                            ▼
                              stdout / --out file (CSV)
```

---

## Key Algorithms

| Component | Algorithm | Reference |
|-----------|-----------|-----------|
| Coordinate frame | WGS-84 ECEF + ENU rotation matrices | ISO 19111 / WGS-84 |
| Track initialization | Weighted Stansfield 3D (WS3D) | Stansfield 1947 / Groves 2013 |
| State estimation | Unscented Kalman Filter (UKF) | Julier & Uhlmann 1997 |
| Process model | Continuous White Noise Acceleration (CWNA) | Bar-Shalom et al. |
| Data association (2-sensor) | Hungarian algorithm | Kuhn-Munkres O(n³) |
| Data association (S>2 sensors) | Lagrangian relaxation MDA | Poore & Rijavec 1993 |
| Gating | Mahalanobis distance + χ²(2 dof) | Bar-Shalom & Li |
| Track management | TOMHT N-scan pruning | Reid 1979 |
| Memory | Zero-allocation ObjectPool | - |

---

## Build

### Dependencies

| Library | Version | Purpose |
|---------|---------|---------|
| Eigen3 | ≥ 3.3 | All linear algebra (header-only) |
| GeographicLib | optional | Hardened WGS-84 transforms |

```bash
# Ubuntu / Debian
sudo apt install libeigen3-dev libgeographic-dev

# macOS
brew install eigen geographiclib
```

### Compile

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### Run tests

```bash
cd build && ctest --output-on-failure
# or directly: ./mtt_test
```

---

## Input Format

Space- or comma-separated columns, one measurement per line:

```
timestamp_s  azimuth_deg  elevation_deg  sensor_lat_deg  sensor_lon_deg  sensor_alt_m
```

Lines starting with `#` are comments.

Example:

```
# timestamp   az     el     lat        lon        alt
0.000000      12.5   8.3    51.477900  -0.001500  46.0
0.000000      14.1   8.7    51.478100   0.012300  46.0
0.033333      12.6   8.4    51.477900  -0.001500  46.0
0.033333      14.2   8.6    51.478100   0.012300  46.0
```

---

## Output Format

CSV with header:

```
timestamp,track_id,state,px_m,py_m,pz_m,vx_mps,vy_mps,vz_mps
```

| Column | Type | Description |
|--------|------|-------------|
| `timestamp` | float | Time of last update [s] |
| `track_id` | uint64 | Unique persistent track identifier |
| `state` | int | 0=Tentative, 1=Confirmed, 2=Coasting |
| `px/py/pz` | float [m] | ECEF position |
| `vx/vy/vz` | float [m/s] | ECEF velocity |

---

## Command-Line Options

```
Usage: mtt_tracker [options] [input_file]

  --gate N     Chi-squared gating threshold (default 13.816 = 99.9% for 2 DoF)
  --psd N      Process noise power spectral density m²/s³ (default 3.0)
  --noise N    Measurement noise variance per axis deg² (default 0.01)
  --confirm N  Consecutive hits required to confirm track (default 3)
  --delete N   Consecutive misses before track deletion (default 5)
  --nscan N    N-scan pruning history depth (default 5)
  --verbose    Print LLA+speed summaries to stderr
  --out FILE   Write output to FILE instead of stdout
```

### Tuning guidance

**Dense bird/aircraft mix** (default):
```bash
./mtt_tracker --psd 3.0 --gate 13.816 --confirm 3 input.txt
```

**Heavy clutter (pollen/insects), clean aircraft only:**
```bash
./mtt_tracker --psd 1.0 --gate 9.21 --confirm 4 input.txt
```

**Highly maneuvering targets:**
```bash
./mtt_tracker --psd 10.0 --gate 13.816 input.txt
```

---

## Module Overview

```
include/mtt/
├── types.hpp            WGS-84 constants, Eigen aliases, MHTConfig
├── geodesy.hpp          LLA↔ECEF, ENU axes (+ src/geodesy.cpp)
├── sensor_manager.hpp   Thread-safe sensor registry
├── observation_model.hpp h(x) bearing function, Jacobian, conversion
├── ukf.hpp              Unscented Kalman Filter (predict + update)
├── ws3d.hpp             Weighted Stansfield 3D triangulator
├── track.hpp            Track object: UKF state, scoring, lifecycle
├── object_pool.hpp      Zero-allocation fixed-capacity pool
├── gating.hpp           Tier-1 Mahalanobis gating
├── assignment.hpp       Hungarian (2D) + Lagrangian MDA (S-D)
├── tomht.hpp            Full TOMHT pipeline
├── blocking_queue.hpp   Thread-safe bounded queue
├── input_parser.hpp     Text-stream parser
└── output_serializer.hpp CSV output formatter
```

---

## Design Decisions vs. Document Spec

| Spec requirement | Implementation |
|-----------------|----------------|
| WGS-84 ECEF state space | ✅ 6D [px,py,pz,vx,vy,vz] in ECEF |
| Degrees → radians conversion with σ² scaling | ✅ `convert_measurement()` |
| UKF over EKF (zenith singularity avoidance) | ✅ Merwe scaled sigma-points |
| TOMHT with persistent track IDs | ✅ Atomic uint64 counter |
| Coarse Mahalanobis gating (Tier 1) | ✅ `CoarseGate` |
| S-D assignment via Lagrangian relaxation (Tier 2) | ✅ `LagrangianMDA` |
| WS3D initialization | ✅ 2D planar + 1D height |
| Asynchronous sequential fusion | ✅ dt-based UKF predict per measurement |
| Zero-allocation ObjectPool | ✅ 10k pre-allocated tracks |
| Multithreaded I/O + tracker + output | ✅ 3 dedicated threads |
| GeographicLib integration | ✅ CMake auto-detect, fallback Bowring |
| Eigen3 with SIMD | ✅ Linked via find_package, -march=native |
| CWNA process model | ✅ F and Q matrices |
