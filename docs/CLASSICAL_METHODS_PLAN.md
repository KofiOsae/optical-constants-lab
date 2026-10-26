# Classical optical methods — implementation specification

## Scope and dependency order
C1 (material CSV) → C2 (validated coherent TMM) → C3 (forward spectral simulator) → C4 (dispersion models) → C5 (selective multilayer inversion) → C6 (guided labs and export).
Do not describe a module as validated until its automated tests have passed.

## C1 — Material data library
- Accept CSV with wavelength_um,n,k and wavelength_nm,n,k; support n-only by explicit k=0 choice, never silently.
- Accept tab-separated files after explicit format detection; reject duplicate, nonfinite, nonmonotonic or nonpositive wavelengths and negative k.
- Preserve source DOI/URL, reference, digitization method, temperature, uncertainty, filename, units, and interpolation choice.
- Use bounded interpolation of n and k on overlapping wavelength grids; do not extrapolate without explicit warning.
- Export normalized wavelength_um,n,k CSV and a provenance metadata JSON.

## C2 — TMM
- Independent forward solver in multilayer.py: arbitrary finite coherent isotropic planar layers with complex wavelength-dependent n+ik, thickness nm, incidence angle and s/p.
- Distinguish lossless incident/exit half-spaces from absorbing finite layers. Explicitly document incoherent substrate/backside reflection as unsupported until modeled.
- Tests: zero-layer Fresnel, one-layer equivalence, zero-thickness, quarter-wave stack, TE/TM oblique incidence, R+T=1 for lossless stacks, R+T<=1 for passive absorbers, convergence/stability of thick absorbers.
- Cross-compare with FMM Research Workbench TMM on identical geometry, optical constants, angle and coherence conventions; save differences and tolerances.

## C3 — Forward simulator UI
- Layer table: order, name, thickness nm, material source (constant, CSV, dispersion model), fixed/fitted toggle.
- Upload material CSV; show wavelength overlap and interpolation preview; export R,T,A and layer stack configuration.
- Simulate 1-3 layer examples first; do not add arbitrary stack fitting until forward tests pass.

## C4 — Dispersion models
- Cauchy and Sellmeier: transparent-region models; avoid using beyond validity and poles.
- Lorentz: epsilon_inf and 1-3 bound-charge oscillators with explicit E0, strength and gamma units.
- Drude-Lorentz: free-carrier plus bound-charge contributions; test causality and passivity conventions.
- Tauc-Lorentz: later, with explicitly documented material/edge assumptions and KK treatment.
- Benchmark parameter recovery using noise-free then noisy synthetic data, multiple initial seeds, noise-weighted residuals, parameter error, n/k RMSE, Jacobian and holdout prediction.

## C5 — Inversion
- Fit selected unknown layers, fixed known materials; regularized spline option and bounded causal parametric models.
- Preserve measurement geometry, noise estimates, instrument R/T conventions and substrate treatment.
- Distinguish training from excluded-film predictions; never tune against holdout.

## C6 — Course and reproducibility
- Tutorial examples: air/glass Fresnel, single film, anti-reflection quarter-wave coating, dielectric Bragg mirror, absorbing organic layer, film on silicon (only after substrate model verified).
- Each example: assumptions, derivation, downloadable source n/k CSV, stack JSON, simulated R/T CSV, expected numerical checks and exercises.
- Report numerical model version, settings, source provenance and uncertainty; literature extraction is not independent experimental verification.

## Operational priorities
- Fix plotting labels; execute tests/test_lorentz_forward.py and tests/test_multilayer.py in CI.
- Profile slow two-stage spline fit; retain numerical objective and settings.
- Free-hosting job design must survive web timeouts; in-memory threads alone are not durable. Identify persistence and restart strategy before claiming background jobs are reliable.
