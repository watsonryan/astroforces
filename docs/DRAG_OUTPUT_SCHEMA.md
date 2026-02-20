# Drag Batch Output Schema

Schema id: `drag_batch_v1`

## Input CSV
Each non-header row:
- `epoch_utc_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps`

## CSV Output
First line is metadata as a comment:
- `#record_type=metadata,...`

Then header row:
- `epoch_utc_s`
- `ax_mps2`, `ay_mps2`, `az_mps2`
- `rho_kg_m3`, `temp_k`
- `vrel_mps`, `q_pa`
- `area_m2`, `cd_eff`
- `wx_source`, `wx_interp`, `wx_extrap`
- `f107`, `f107a`, `ap_daily`, `kp_daily`, `ap_3h`, `kp_3h`
- `status`

## JSON Output (NDJSON)
Line 1:
- metadata object:
  - `record_type: "metadata"`
  - `schema: "drag_batch_v1"`
  - `project`, `generated_unix_utc`, `model`, `model_data`, `wind`, `wind_data`, `weather_csv`

Subsequent lines:
- sample objects:
  - `record_type: "sample"`
  - `schema: "drag_batch_v1"`
  - same physical fields as CSV output
