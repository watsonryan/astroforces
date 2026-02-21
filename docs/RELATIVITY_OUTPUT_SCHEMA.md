# Relativity Output Schema

## Single-State CLI
Command:
- `relativity_cli x_eci_m y_eci_m z_eci_m vx_eci_mps vy_eci_mps vz_eci_mps epoch_utc_s jpl_ephemeris_file [use_geodesic]`

Output fields:
- Acceleration: `ax`, `ay`, `az`, `amag`
- Term magnitudes:
  - `schw` (spherical central-body term)
  - `geo` (geodesic precession)
  - `lt` (Lense-Thirring)
  - `j2` (relativistic oblateness/J2 correction)
  - `rot` (rotational-energy correction)

## Notes
- The implementation follows Eq. 25 term decomposition with configurable PPN parameters.
- The rotational-energy term is implemented as a separate switch, though in practice it is strongly coupled with the recovered effective `J2`.

