# V2 guided example library — synthetic, independently checkable inputs

**Provenance:** The following datasets are generated analytically using the exact normal-incidence coherent Air / thin film / semi-infinite transparent substrate Fresnel slab formula. No laboratory measurements or literature-derived optical constants are implied. Uniform wavelength sampling: 450–850 nm, 5 nm steps; reflectance and transmittance are absolute fractions, not percentages. Substrate n=1.46, k=0. Ambient n=1. The files have no added measurement noise.

**Test A — Transparent single film (120 nm):** [CSV](transparent_film_120nm.csv). Film n=1.65 constant and k=0, thickness 120 nm. Guide: choose Single Spectrum → upload → preview → Cauchy (transparent) if available, or Sellmeier with appropriate fixed poles. Set substrate n=1.46 and film thickness 120 nm; disable mixing and noise. Check R+T≈1. A single finite-band spectrum can have ambiguous index/thickness solutions: don't claim unique recovery from optimizer convergence.

**Test B — Transparent joint films (120 nm + 240 nm):** [120 nm](transparent_film_120nm.csv), [240 nm](transparent_film_240nm.csv). Both have n=1.65, k=0; different known thicknesses. Guide: choose Joint Spectrum → upload BOTH → set thicknesses explicitly 120 and 240 nm → Cauchy/Sellmeier transparent model → hold out none → fit. Compare common reconstructed n,k and R/T per film. This is an identifiability demonstration; inspect residuals and recovered parameters, not RMSE alone.

**Test C — Spectrally absorbing film (180 nm):** [CSV](absorbing_film_180nm.csv). Film thickness 180 nm. Synthetic optical constants: n(λ)=1.73+0.12(550/λ−1); k(λ)=0.07 exp[−0.5((λ−610)/60)^2], λ in nm. Guide: choose Single Spectrum with a flexible dispersion model or a suitable causal representation, set thickness 180 nm, substrate 1.46, no mixing/no noise initially. Since this demonstration uses a descriptive n,k function rather than a globally causal material model, **do not treat it as a KK/SSKK causality benchmark**. Evaluate R/T reconstruction and fit residuals. A particular causal fit cannot be guaranteed to exactly recover this deliberately noncausal toy dispersion.

## Model-selection flow, practical checklist
1. Define geometry: air/film/substrate, wavelength units, substrate n,k, known thicknesses.
2. Confirm spectral content and physical ranges in preview, R≥0, T≥0, passive energy balance.
3. Transparent sample → Cauchy/Sellmeier. Absorbing film → causal oscillator fit if physically motivated, otherwise flexible spline as exploratory fit.
4. Start with fixed thickness, zero mixing and zero noise. Fit an informative wavelength band, then expand.
5. Inspect R/T overlay, spectral residuals, recovered n,k, parameter bounds, and independent references where available. For causal models inspect finite-band KK/SSKK with the documented anchor assumptions.
6. Export project and methods package. Repeat from scratch and compare; note solver/model assumptions.

## Required numerical verification
Run a separate forward Fresnel/TMM implementation to compare each CSV; report max |ΔR| and |ΔT|. Then perform inversion and report both fit residual and deviation from the known synthetic n,k. Do not claim that these demonstration files represent measured data. Report multiple initializations and perturbation tests.

For v2, these are **teaching/example datasets**. Independent published dataset onboarding and a held-out benchmark suite remain further work.
