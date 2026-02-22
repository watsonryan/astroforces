# Relativity Model Notes (Eq. 25)

This implementation follows the Eq. 25 decomposition from:
- Decaia et al., Section 6.8, Eq. 25 decomposition:
  `https://www.iag-ggos.org/sci/itrsstandard_docs/decaiaacos.pdf`
  (ITRS standards technical note document used for force-model term decomposition).

Implemented additive terms:
- Spherical central-body term (PPN `beta`, `gamma`)
- Geodesic precession term (de Sitter)
- Lense-Thirring term (frame dragging)
- Oblateness/J2 relativistic correction
- Rotational-energy correction

Configuration knobs:
- `ppn_beta`, `ppn_gamma`, `lense_thirring_parameter`
- Earth/Sun `mu`, `c`, `J2`, Earth reference radius
- Earth spin unit vector, angular momentum per unit mass (`Je`), rotational energy per unit mass (`Te`)

Notes:
- Geodesic precession requires ephemerides (`Earth wrt Sun` state).
- The rotational-energy term is modeled separately, though standards literature notes practical coupling with recovered effective `J2`.
- For default GR behavior: `beta = gamma = 1`, `lense_thirring_parameter = 1`.
