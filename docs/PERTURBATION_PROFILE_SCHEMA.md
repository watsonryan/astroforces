# Perturbation Profile Output Schema

Schema id: `perturbation_profile_v1`

## CLI
`perturbation_profile_cli` arguments:
- `output_csv`
- `alt_min_km`
- `alt_max_km`
- `samples`
- `dtm_coeff_file`
- `space_weather_csv`
- `jpl_ephemeris_file` (optional)
- `epoch_utc_s` (optional)

## CSV Output
Header starts with:
- `altitude_km`

Then one or more component columns:
- `<component_name>_mps2`

Current default components:
- `drag_mps2`
- `third_body_sun_mps2` (if ephemeris provided)
- `third_body_moon_mps2` (if ephemeris provided)

Final columns:
- `total_mps2` (root-sum-square of component magnitudes)
- `status` (`astroforces::atmo::Status` enum as integer)

## Extensibility
New perturbation components can be added without schema break by appending additional
`<component_name>_mps2` columns before `total_mps2,status`.

