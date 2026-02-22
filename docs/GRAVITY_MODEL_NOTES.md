# Gravity SPH Model Notes

## Scope

`GravitySphAccelerationModel` evaluates Earth gravity as:
- degree-0 central term
- non-central spherical harmonics from `Cnm/Snm` coefficients (`n>=2`) to configurable degree/order
- tide stack terms compatible with the Ginan-style decomposition:
- solid Earth tide-1 (Sun/Moon body tide)
- solid Earth tide-2 (frequency-dependent `C20/C21/S21/C22/S22`)
- pole tide (solid)
- pole tide (ocean)
- constituent ocean tide
- constituent atmospheric tide

The model accepts state in `ECI` or `ECEF`.
For `ECI`, a simple GMST z-rotation is used internally to evaluate the Earth-fixed field.

## Coefficient File

The gravity model reads ICGEM-like `.gfc` data:
- header fields: `earth_gravity_constant`, `radius`, `tide_system`
- data rows: `gfc n m Cnm Snm ...`

If the source is not tide-free and `convert_to_tide_free=true`, `C20` is corrected to tide-free compatibility.

Important for EIGEN-6S4 (Version 2):
- EIGEN files may contain time-variable records (`gfct`, `trnd`, `acos`, `asin`).
- This implementation synthesizes these time-variable terms when present.

Reference source for model download/comparison:
- `https://icgem.gfz.de/tom_longtime`

## Tide Inputs

- Sun/Moon solid tides require JPL ephemeris.
- Pole tides require IERS finals EOP (`xp/yp`).
- Ocean pole tide can optionally use a coefficient file (`n m cnmp cnmm snmp snmm`);
  otherwise the built-in IERS-style fallback coefficients are used.
- Constituent ocean/atmospheric tides require potential coefficient files in the expected constituent format.

## SPH Acceleration Synthesis

The non-central acceleration follows the same derivative synthesis pattern used in Ginan (`accelSPH`):
- build normalized associated Legendre functions `Pnm` and derivatives `dPnm`
- accumulate radial/latitudinal/longitudinal potential derivatives
- scale by `mu` and project to Cartesian acceleration

## Configuration Notes

- `max_degree` defaults to `360` and is clamped to the coefficient file maximum.
- no numerical Jacobians are used.
- outputs are in m/s^2 and returned by term (`central`, `sph`, `solid_tide_sun`, `solid_tide_moon`, `solid_tide_freqdep`, `pole_tide_solid`, `pole_tide_ocean`, `ocean_tide`, `atmos_tide`).
