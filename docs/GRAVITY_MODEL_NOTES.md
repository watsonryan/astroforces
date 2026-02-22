# Gravity SPH Model Notes

## Scope

`GravitySphAccelerationModel` evaluates Earth gravity as:
- degree-0 central term
- non-central spherical harmonics from `Cnm/Snm` coefficients (`n>=2`) to configurable degree/order
- optional solid Earth tide correction (IERS2010 Eq. 6.6/6.7 style coefficient deltas from Sun/Moon)

The model accepts state in `ECI` or `ECEF`.
For `ECI`, a simple GMST z-rotation is used internally to evaluate the Earth-fixed field.

## Coefficient File

The gravity model reads ICGEM-like `.gfc` data:
- header fields: `earth_gravity_constant`, `radius`, `tide_system`
- data rows: `gfc n m Cnm Snm ...`

If the source is not tide-free and `convert_to_tide_free=true`, `C20` is corrected to tide-free compatibility.

Important for EIGEN-6S4 (Version 2):
- EIGEN files may contain time-variable records (`gfct`, `trnd`, `acos`, `asin`).
- Current implementation uses the static `gfc` coefficients only.
- Time-variable gravity terms are not yet synthesized in this repo.

Reference source for model download/comparison:
- `https://icgem.gfz.de/tom_longtime`

## SPH Acceleration Synthesis

The non-central acceleration follows the same derivative synthesis pattern used in Ginan (`accelSPH`):
- build normalized associated Legendre functions `Pnm` and derivatives `dPnm`
- accumulate radial/latitudinal/longitudinal potential derivatives
- scale by `mu` and project to Cartesian acceleration

## Solid Earth Tides

When enabled:
- Sun/Moon positions are obtained from `jplEphem`
- each body contributes `dCnm/dSnm` corrections (degrees 2-4) using elastic Love numbers
- corrected coefficients are then synthesized by the same SPH evaluator

## Configuration Notes

- `max_degree` defaults to `360` and is clamped to the coefficient file maximum.
- no numerical Jacobians are used.
- outputs are in m/s^2 and returned by term (`central`, `sph`, `solid_tide_sun`, `solid_tide_moon`).
