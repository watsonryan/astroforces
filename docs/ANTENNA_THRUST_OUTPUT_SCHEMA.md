# Antenna Thrust Output Schema

Schema id: `antenna_thrust_v1`

## CLI
`antenna_thrust_cli` arguments:
- `x_eci_m`
- `y_eci_m`
- `z_eci_m`
- `vx_eci_mps`
- `vy_eci_mps`
- `vz_eci_mps`
- `epoch_utc_s`
- `mass_kg` (optional)
- `transmit_power_w` (optional)
- `efficiency` (optional)
- `mode` (optional): `velocity|nadir|custom_eci|body_fixed`
- `dir_x` (optional; used by `custom_eci` or `body_fixed`)
- `dir_y` (optional; used by `custom_eci` or `body_fixed`)
- `dir_z` (optional; used by `custom_eci` or `body_fixed`)

## Batch CLI
`antenna_thrust_batch_cli` arguments:
- `input_csv`
- `output_csv`
- `mass_kg` (optional)
- `transmit_power_w` (optional)
- `efficiency` (optional)
- `mode` (optional): `velocity|nadir|custom_eci|body_fixed`
- `dir_x` (optional)
- `dir_y` (optional)
- `dir_z` (optional)

Input row format:
- `epoch_utc_s,x_eci_m,y_eci_m,z_eci_m,vx_eci_mps,vy_eci_mps,vz_eci_mps`

## Batch CSV Output
Header:
- `epoch_utc_s`
- `ax_mps2`
- `ay_mps2`
- `az_mps2`
- `amag_mps2`
- `thrust_n`
- `effective_power_w`
- `mass_kg`
- `dir_x`
- `dir_y`
- `dir_z`
- `status`

Status is `astroforces::core::Status` enum integer.

