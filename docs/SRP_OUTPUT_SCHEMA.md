# SRP Output Schema

## Single-State CLI
Command:
- `srp_cli x_eci_m y_eci_m z_eci_m vx_eci_mps vy_eci_mps vz_eci_mps epoch_utc_s jpl_ephemeris_file [mass_kg] [area_m2] [cr] [use_eclipse]`

Output fields:
- Acceleration: `ax`, `ay`, `az`, `amag`
- Diagnostics: `p_pa`, `r_sun_m`, `area`, `cr`, `eclipsed`

## Batch CLI
Command:
- `srp_batch_cli input_csv output_csv jpl_ephemeris_file [mass_kg] [area_m2] [cr] [use_eclipse]`

Input row format:
- `epoch_utc_s,x_eci_m,y_eci_m,z_eci_m,vx_eci_mps,vy_eci_mps,vz_eci_mps`

Output CSV columns:
- `epoch_utc_s`
- `ax_mps2`, `ay_mps2`, `az_mps2`, `amag_mps2`
- `solar_pressure_pa`, `sun_distance_m`
- `area_m2`, `cr`, `eclipsed`
- `status` (`astroforces::core::Status` enum value as integer)

