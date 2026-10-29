# Organic dye dispersion models — implementation and release checklist

This branch begins a **shared** dispersion-engine extension. It is NOT a completed three-workflow integration or a release candidate.

## Model conventions
* Wavelength input: nm; photon energy: 1239.841984 / wavelength_nm (eV).
* Dielectric convention: exp(-i omega t), with N = sqrt(epsilon) = n + i k.
* Positive oscillator strengths and damping; k >= 0 for passive materials.
* The implemented Gaussian-distributed Lorentz representation is causal because each component is a passive Lorentz pole. It is not an arbitrary Gaussian applied to k.
* Modified Lorentz currently means a positive narrow/broad Lorentz mixture, not a unique literature convention.
* Vibronic progression uses six truncated Poisson-weighted sidebands; this phenomenological choice requires explicit documentation and testing.

## Shared fit specification (must be completed before release)
1. Centralize model selection, parameter layouts, bounds, serialization, and warm-start compatibility.
2. Wire all applicable models to **single-spectrum, joint-spectrum and multilayer** inversions using one canonical dielectric function and compatible result schemas.
3. Retain phase/envelope, Cauchy and Sellmeier as transparent-region-only specializations; do not treat them as general absorbing dye fits.
4. Implement UI selectors and controls for all three workflows. Ensure exact same model names, oscillator limits, default bounds, wavelength-range fitting, noise support, plots, reference overlays, and model exports. Disable features with clear explanations only where physically inapplicable.
5. Keep training vs held-out metrics separate; record model, oscillator count, coefficients, bounds, fit domain, thickness choices, uncertainties and termination.
6. Compare with Lorentz baseline and synthetic ground truth, including R/T residual shape and accuracy of recovered n and k; test multiple film thicknesses and a held-out sample.
7. Benchmark convolution quadrature runtime and resolution; optimize before deployment to free hosting.
8. Update TUTORIAL.md and the embedded tutorial in web/workspace.html, plus multilayer help and README. Include a reproducible worked organic-dye example and limitations.
9. Check workbook export: no claims that novel models reconstruct in Excel until model-specific formulas and tests exist.
10. Run CI and perform an end-to-end deployed preview, fit and export across all three workflows before tagging a release.

## Status
* Shared organic dielectric evaluator: implemented in organic_dispersion.py.
* Unit tests: added, not yet confirmed by CI.
* Optimizer integration, UI parity, tutorial updates and deployment: pending. 
