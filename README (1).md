# Optical Constants Inversion: standalone reference engine

**Status: NOT integrated into FMM Research Workbench or pushed to GitHub.**

## Supported now
- Coherent, uniform single film on a semi-infinite, nonabsorbing substrate.
- Air or other transparent incident medium; s/p polarization and oblique incidence.
- Joint fitting of multiple R/T spectra sharing a common wavelength-dependent `n,k` spline.
- Fixed, supplied film thickness per spectrum, and common `n,k` across samples.
- Sensitivity/Jacobian and noise-weighted Fisher information.
- Basic physics and synthetic recovery tests.

## Not yet supported (required for full Workbench integration)
- Original substrate-dependent Swanepoel envelope implementation; phase-I implementation should be ported and benchmarked separately.
- Coherent-incoherent finite glass substrates, multilayer TMM, roughness/scattering.
- Fit of unknown film thickness, KK-consistent models, ellipsometry, instrumental correlation.
- Data import UI, CSV validation, uncertainty bootstrap, automatic model selection, Render integration.

## Multi-spectrum measurement metadata
Each spectrum requires wavelength in **nm**, R/T as **power fractions**, known thickness in **nm**, substrate index, incidence angle in degrees, polarization, and uncertainty standard deviations. All spectra must correspond to **the same material optical constants**; different concentrations, deposition conditions, porosity, or film morphology violate the shared-material assumption.

Spectra may have different wavelength grids. If wavelengths do not overlap, joint fitting cannot by itself guarantee trustworthy interpolation across gaps. Do not extrapolate the spline outside the collective fit domain.

## Installation and tests
`pip install numpy scipy pytest`
`pytest -q tests/`

## Example
```python
import numpy as np
from inversion import Spectrum, coherent_film, fit_joint
wl=np.linspace(400,800,81)
samples=[]
for d,ns in [(50,1.0),(1000,1.0),(50,1.5),(1000,1.5)]:
    R,T=coherent_film(wl,np.full(len(wl),1.8+.025j),d,ns)
    samples.append(Spectrum(wl,R,T,d,ns))
fit, optical_constants=fit_joint(samples,n_basis=10)
n,k=optical_constants(wl)
print(fit.success, np.mean(n), np.mean(k))
```

## Workbench integration checklist
1. Inspect existing routes, Material database and Fresnel/TMM interfaces.
2. Wrap engine as a backend service with explicit units and configuration schema.
3. Add data import, validation and multi-spectrum sample table to UI.
4. Add original Swanepoel method and distinct spline/Fresnel and causal model options.
5. Plot per-spectrum measured vs fitted R/T, residuals, recovered n/k and uncertainties.
6. Add physics limits, ground-truth recovery, noise, holdout, and regression tests.
7. Run against actual Render environment and verify server-side limits.

This is a reference implementation, not a finished production app.
