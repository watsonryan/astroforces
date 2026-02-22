# Perturbation Profile Output Schema

Schema id: `perturbation_profile_v2`

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
- `gravity_gfc_file` (optional)
- `gravity_max_degree` (optional)
- `eop_finals_file` (optional)
- `ocean_tide_file` (optional)
- `atmos_tide_file` (optional)
- `ocean_pole_tide_file` (optional)

## CSV Output
Header starts with:
- `altitude_km`

Then one or more component columns:
- `<component_name>_mps2`

Current default components:
- `drag_mps2`
- `antenna_thrust_mps2`
- `gravity_central_mps2` (if gravity file provided)
- `gravity_sph_tides_mps2` (if gravity file provided)
- `gravity_tide_solid_sun_mps2` (if gravity file + ephemeris provided)
- `gravity_tide_solid_moon_mps2` (if gravity file + ephemeris provided)
- `gravity_tide_solid_freqdep_mps2` (if gravity file provided)
- `gravity_tide_pole_solid_mps2` (if gravity file + EOP provided)
- `gravity_tide_pole_ocean_mps2` (if gravity file + EOP provided)
- `gravity_tide_ocean_mps2` (if gravity file + ocean tide file provided)
- `gravity_tide_atmos_mps2` (if gravity file + atmos tide file provided)
- `earth_radiation_mps2`
- `relativity_mps2`
- `srp_mps2` (if ephemeris provided)
- `third_body_sun_mps2` (if ephemeris provided)
- `third_body_moon_mps2` (if ephemeris provided)

Final columns:
- `total_mps2` (root-sum-square of component magnitudes)
- `status` (`astroforces::core::Status` enum as integer)

## Extensibility
New perturbation components can be added without schema break by appending additional
`<component_name>_mps2` columns before `total_mps2,status`.
