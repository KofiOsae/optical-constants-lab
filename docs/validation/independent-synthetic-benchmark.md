# Independent synthetic benchmark — Optical Constants Laboratory v2

## Scope and provenance

The four CSV files in `examples/v2` are synthetic **teaching** spectra (450–850 nm, 5 nm steps). They are not measured R/T, published material constants, or independent evidence of accurate optical-constant recovery. `tests/test_independent_fresnel_reference.py` is an **independent implementation** of normal-incidence amplitude recursion; it does not import `multilayer_forward` or the inversion engine. It evaluates forward calculations from explicitly known film indices and thicknesses. It does **not** verify parameter identifiability, global causality or optimizer convergence.

Each file contains columns `wavelength_nm,R,T`, with `R,T` as absolute fractions. Assumed ambient index 1.0, nonabsorbing semi-infinite substrate 1.46, normal incidence, coherent films, and no instrument noise.

| Test | Stack, incident to substrate | Expected optical constants |
| --- | --- | --- |
| Single transparent | 120 nm transparent film | n = 1.65, k = 0 |
| Joint transparent | Two **separate samples**, 120 and 240 nm of the same film | n = 1.65, k = 0 |
| Single absorbing | 180 nm film | n(lambda) = 1.73 + 0.12(550/lambda - 1); k(lambda) = 0.07 exp[-0.5((lambda - 610)/60)^2], lambda in nm |
| Multilayer cap | 50 nm known cap followed by 150 nm unknown transparent film | cap n = 1.80, k = 0; film n = 1.65, k = 0 |

The absorbing example uses a descriptive finite-band function, **not** an independently causal dielectric response; do not use it to grade KK/SSKK reconstruction.

## Reference equations

Let `N_j = n_j + i k_j` use the solver's absorption convention. For light incident from medium `j` to `j+1` at normal incidence,

```text
r_j = (N_j - N_(j+1))/(N_j + N_(j+1))
t_j = 2 N_j/(N_j + N_(j+1))
p = exp(i 2 pi N_(j+1) d_(j+1)/lambda)
r_eff = (r_j + r_next p^2)/(1 + r_j r_next p^2)
t_eff = (t_j t_next p)/(1 + r_j r_next p^2)
R = |r_eff|^2
T = Re(N_substrate)/Re(N_ambient) |t_eff|^2
```

The recursion is evaluated from the substrate back toward the incident medium. The test computes the reflectance and transmittance with independent equations and checks against the CSV values, passivity, and `R+T=1` for transparent samples.

## Verification criteria

Run `python -m pytest -q tests/test_independent_fresnel_reference.py`. For all four files require: strictly increasing wavelengths, 81 rows, max absolute difference `|Delta R| < 1e-8` and `|Delta T| < 1e-8`, R,T nonnegative, passive `R+T<=1+1e-8`; for transparent samples require `max|R+T-1|<1e-8`. For the absorbing spectrum, require positive absorption somewhere.

An independently evaluated complex-amplitude Fresnel recursion (implemented separately from the laboratory code) produced maximum absolute differences below `5e-10` for each of the three tested example types: transparent 120 nm, absorbing 180 nm, and multilayer. The repository test also covers transparent 240 nm; the CI job is the authority for its actual result.

## Inverse recovery acceptance: further validation required

To make this a genuine **inverse** benchmark, start from blind (non-truth) initial values, preserve known thickness for the basic case, and use repeated optimizer starts. Export:
1. recovered `n(lambda)`, `k(lambda)` versus known truth; report `MAE(n)`, `MAE(k)`, maximum absolute error and systematic bias over the measured wavelengths;
2. `R/T` residuals for each sample, including spectral plots and weighted RMSE with stated weighting convention;
3. parameter recovery errors, any active bounds, optimizer termination, rank/condition number where meaningful, and variation across starts;
4. perturbations of thickness and substrate index; noise trials with fixed seed and at least three signal/noise levels;
5. independent held-out thickness/film spectrum omitted from *both* optimizer fitting and hyperparameter selection.

Report a failed or ambiguous recovery as such, even when spectral RMSE is tiny. Do not assign numerical pass thresholds before reference recovery distributions and measurement uncertainties have been determined.

## Reproduction information

Record the exact application commit, Python version, NumPy/SciPy versions, input-file hashes, model family, chosen settings and bounds, random seed, fit wavelength window, convergence diagnostics and archived CSV outputs. Never claim a synthetic demo has established experimental accuracy.
