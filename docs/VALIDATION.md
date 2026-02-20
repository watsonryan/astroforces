# Validation Plan

## Stage 1 (Current)
- Analytic drag-equation unit tests
- End-to-end smoke test with deterministic exponential atmosphere
- Macro-model attitude projection tests for area and effective Cd
- Space-weather parser regression tests (synthetic and real-snippet CelesTrak rows)

## Stage 2 (Integration)
- Adapter parity tests versus source model repos
- Cross-model consistency checks on shared scenarios
- NRLMSIS Ap-history adapter plumbing validation
- Drag batch output schema checks (CSV/JSON rows)
- Direct adapter-vs-model parity assertion for NRLMSIS `ap_history` mode

## Stage 3 (Operational)
- Golden trajectory drag acceleration regressions
- Performance budget checks in CI
