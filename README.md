# codex-test

MTT is a C++ scaffold for real-time multi-target tracking from azimuth/elevation measurements produced by 2-4 ground cameras.

## Build
```bash
cmake -S . -B build
cmake --build build -j
```

## Run
```bash
./build/mtt_cli measurements.csv
```

Input format (CSV):

```text
timestamp, azimuth, elevation, sensor latitude, sensor longitude, sensor altitude
```

The tracker estimates and prints one line per track with:
- `track_id`
- ECEF position `[px, py, pz]` in meters
- ECEF velocity `[vx, vy, vz]` in m/s

See detailed design notes in [`docs/MTT_DESIGN.md`](docs/MTT_DESIGN.md).
