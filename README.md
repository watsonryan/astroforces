# astrodynamics-forces-cpp

Unified C++20 astrodynamics perturbation modeling platform, with drag implemented first.

## Architecture
```mermaid
flowchart LR
  A[atmo-core<br/>types, frames, interfaces] --> B[space-weather<br/>CelesTrak CSV]
  A --> C[sc-props<br/>cannonball + macro surfaces]
  A --> D[adapters<br/>NRLMSIS / DTM2020 / HWM14]
  A --> E[drag-core<br/>DragAccelerationModel]
  A --> F[forces<br/>IPerturbationModel + PerturbationStack]
  B --> E
  C --> E
  D --> E
  E --> F
  T[third-body<br/>Sun+Moon via jplEphem] --> F
  F --> G[apps/drag-cli]
  F --> H[apps/drag_batch_cli]
  F --> P[apps/perturbation_profile_cli]
  P --> Q[scripts/plot_perturbation_profile.py]
  I[external repos via CPM] --> D
  I --> T
  I --> B
```

## Scope
- Shared force interfaces (`IPerturbationModel`, `PerturbationStack`)
- Drag acceleration core (relative velocity + ballistic term)
- Spacecraft geometry/surface property hooks
- Adapter integration points for:
  - NRLMSIS 2.1 (`nrlmsis-2_1`)
  - DTM2020 (`dtm2020`)
  - HWM14 (`hwm14`)

## Build
```bash
cmake --preset macos-debug
cmake --build --preset macos-debug
ctest --preset macos-debug --output-on-failure
```

External model repos are pulled via CPM with HTTPS URLs by default; override cache vars
`DRAGCPP_NRLMSIS21_REPO`, `DRAGCPP_DTM2020_REPO`, and `DRAGCPP_HWM14_REPO` if you prefer SSH.

## CLI
```bash
./build/macos-debug/drag_cli 6778137 0 0 0 7670 0 1000000000
```

Batch CLI:
```bash
./build/macos-debug/drag_batch_cli input.csv output.csv csv nrlmsis /path/to/msis21.parm zero "" /path/to/SW-Last5Years.csv
```
Input row format:
- `epoch_utc_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps`
- Output schema reference: `docs/DRAG_OUTPUT_SCHEMA.md`

Optional weather input:
- Pass CelesTrak 5-year CSV as final arg:
```bash
./build/macos-debug/drag_cli 6778137 0 0 0 7670 0 1000000000 nrlmsis /path/to/msis21.parm zero "" /path/to/SW-Last5Years.csv
```

Weather mapping notes:
- CelesTrak KP columns are parsed in tenths and converted to 0-9 scale.
- `kp_3h_current` is used for DTM delayed 3-hour Kp input.
- Daily `kp` is used for DTM 24-hour mean Kp input.
- `ap_msis_history` is computed from 3-hour AP slots and passed through the NRLMSIS adapter.

Drag area modes:
- Cannonball: set `use_surface_model=false`; area is fixed at `reference_area_m2`.
- Macro-model: set `use_surface_model=true` with surfaces; projected area is computed from plate normals and flow direction in body frame.
- Body-frame flow direction uses `StateVector::body_from_frame_dcm` (row-major DCM).

General perturbation interface:
- Use `astroforces::forces::IPerturbationModel` for each force source.
- Combine models with `astroforces::forces::PerturbationStack`.
- Drag is exposed as `astroforces::drag::DragPerturbationModel` and plugs directly into the same stack used for future gravity/SRP/third-body models.
- Third-body is exposed as `astroforces::forces::ThirdBodyPerturbationModel` (Sun/Moon direct + indirect terms via JPL ephemerides).

Perturbation-vs-altitude profiling:
```bash
./build/macos-debug/perturbation_profile_cli perturbation_profile.csv 200 20000 500 \
  /path/to/operational_regression_coeff.dat \
  /path/to/SW-Last5Years.csv \
  /path/to/linux_p1550p2650.440 \
  1000000000
python3 scripts/plot_perturbation_profile.py perturbation_profile.csv --output-stem perturbation_vs_altitude --column single
```
This generates publication-ready IEEE-style PDF/PNG plots of acceleration magnitude by perturbation component versus altitude.
Output columns include:
- `drag_mps2`
- `third_body_sun_mps2`
- `third_body_moon_mps2`
- `total_mps2`

Notes:
- Third-body columns are component-based and extensible for future bodies (additional `*_mps2` columns).
- Output schema reference: `docs/PERTURBATION_PROFILE_SCHEMA.md`.

Performance benchmark:
```bash
./build/macos-debug/dragcpp_perf_benchmark
```
Environment controls:
- `ASTRO_FORCES_PERF_SAMPLES` (default `40`)
- `ASTRO_FORCES_PERF_ITERS` (default `5000`)
