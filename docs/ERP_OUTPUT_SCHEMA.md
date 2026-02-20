# ERP Output Schema

## Single-State CLI
Command:
- `erp_cli x_frame_m y_frame_m z_frame_m vx_frame_mps vy_frame_mps vz_frame_mps epoch_utc_s [mass_kg] [area_m2] [cr]`

Output fields:
- Acceleration: `ax`, `ay`, `az`, `amag`
- Diagnostics: `p_pa`, `r_earth_m`, `area`, `cr`

## Batch CLI
Command:
- `erp_batch_cli input_csv output_csv [mass_kg] [area_m2] [cr]`

Input row format:
- `epoch_utc_s,x_eci_m,y_eci_m,z_eci_m,vx_eci_mps,vy_eci_mps,vz_eci_mps`

Output CSV columns:
- `epoch_utc_s`
- `ax_mps2`, `ay_mps2`, `az_mps2`, `amag_mps2`
- `earth_radiation_pressure_pa`, `earth_distance_m`
- `area_m2`, `cr`
- `status` (`astroforces::atmo::Status` enum value as integer)

