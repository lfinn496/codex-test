# codex-test

MTT is a C++ scaffold for a real-time multi-target tracker operating on azimuth/elevation measurements from 2-4 ground cameras.

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

See detailed design notes in [`docs/MTT_DESIGN.md`](docs/MTT_DESIGN.md).
