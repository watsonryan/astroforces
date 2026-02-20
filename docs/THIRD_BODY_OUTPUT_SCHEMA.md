# Third-Body Output Schema

## Single-State CLI
Command:
- `third_body_cli x_eci_m y_eci_m z_eci_m vx_eci_mps vy_eci_mps vz_eci_mps epoch_utc_s jpl_ephemeris_file [use_sun] [use_moon]`

Output fields:
- Sun contribution: `sun_ax`, `sun_ay`, `sun_az`, `sun_mag`
- Moon contribution: `moon_ax`, `moon_ay`, `moon_az`, `moon_mag`
- Total contribution: `total_ax`, `total_ay`, `total_az`, `total_mag`

## Batch CLI
Command:
- `third_body_batch_cli input_csv output_csv jpl_ephemeris_file [use_sun] [use_moon]`

Input row format:
- `epoch_utc_s,x_eci_m,y_eci_m,z_eci_m,vx_eci_mps,vy_eci_mps,vz_eci_mps`

Output CSV columns:
- `epoch_utc_s`
- `sun_ax_mps2`, `sun_ay_mps2`, `sun_az_mps2`, `sun_mag_mps2`
- `moon_ax_mps2`, `moon_ay_mps2`, `moon_az_mps2`, `moon_mag_mps2`
- `total_ax_mps2`, `total_ay_mps2`, `total_az_mps2`, `total_mag_mps2`
- `status` (`astroforces::atmo::Status` enum value as integer)

