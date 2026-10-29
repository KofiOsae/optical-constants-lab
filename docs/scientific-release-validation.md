# Scientific release qualification — optical-constant reconstruction

Status: **in progress; not a general accuracy certification**. Owner: Optical Constants Laboratory. Date: 2026-10-29.

## 1. Non-negotiable separation
A low R/T fit error, optimizer success, optical-constant accuracy, identifiability, and out-of-sample predictivity are **five different** results. Never report one as proof of another. Ground truth n,k and held-out measurements must not inform seeds, bounds, regularization, model selection, branch choice, or stopping rules. Pre-register those before the blind inversion. Keep reference data accessible only to the benchmark evaluator.

## 2. Independent evidence tiers
- **Tier A, software parity:** forward-generated observations from the same production model. Necessary for regression, *insufficient* for physical validation.
- **Tier B, independently coded analytic model:** Fresnel/Airy amplitude formulation, with fixed-seed noise. Tests independent implementation and measured-versus-fit consistency within coherent isotropic assumptions. Test: `test_independent_airy_blind_recovery_and_unseen_thickness`.
- **Tier C, cross-solver synthetic:** spectra exported from FMM/RCWA using wavelength-dependent literature RI, with input material RI files and geometry stored separately. Includes Air/Au(30 nm)/fused silica, Air/R6G(80 nm)/Au(30 nm)/silica and Air/R6G(160 nm)/Au(30 nm)/silica. No reference RI is supplied to reconstruction. Must compare the forward parity of *known* layers before attributing errors to the inverse solver.
- **Tier D, experimental:** calibrated measured R/T at known geometry, with independently characterized thickness and trusted RI or ellipsometry comparator, including uncertainty and provenance.

## 3. Benchmark matrix
| ID | Specimen | Model family | Training data | Held out | Intended failure mode |
|---|---|---|---|---|---|
| Q01 | Weakly absorbing transparent film | B-spline, Phase I/II | Two thicknesses | Third thickness | Fringe aliasing, thickness/index coupling |
| Q02 | Single Lorentz resonance | Lorentz and spline | Two thicknesses | Third thickness, angle | Missing resonance and overfitting |
| Q03 | Two Lorentz resonances | One- and two-pole Lorentz | Several spectra | Unseen polarization/angle | Under-modeling and extra degrees of freedom |
| Q04 | Drude-like metal | Drude–Lorentz, B-spline | Au film spectra | Second thickness | Incorrect metal branch and sensitivity |
| Q05 | R6G on Au | B-spline, Lorentz | 80 and 160 nm R6G / known Au | Third R6G thickness or independent angle | Ambiguity and high absorption |
| Q06 | Transparent lossless stack | Any supported stack model | Normal incidence | TE/TM at oblique angle | Incorrect layer ordering and polarization |
| Q07 | Mis-specified Au RI, roughness, incoherent substrate | Relevant fits | R/T | Independent R/T | False numerical success under forward-model mismatch |
| Q08 | Missing R or T, sparse wavelengths, corrupted units | Applicable models | Altered input | None | Graceful rejection or explicit limitation |
| Q09 | Multiple nearly equivalent branches | B-spline | One thin film | Additional thickness | Degeneracy detection without ground truth |
| Q10 | Substrate-index/thickness perturbations | Applicable models | Joint measurements | Independent stack | Bias/uncertainty sensitivity |

Each test must report the exact physical configuration, wavelength grid, model, measurement uncertainty, noise seed, bounds, starting points, optimizer budget, and runtime.

## 4. Pre-registered evaluation (do not tune on held-out data)
For measured samples S, define per-sample spectral RMSE as `sqrt(mean(concatenate([R_pred-R_obs, T_pred-T_obs])**2))` using only available channels. Use a separately labeled weighted RMSE from the actual inverse residual vector. Compute n_RMSE and k_RMSE on the *measured* wavelength interval against independent reference. Also calculate n/k MAE, 95th-percentile absolute error, minimum/maximum and fraction of optical samples hitting bounds, R/T holdout RMSE, and energy-conservation/passivity violations.

Initial *research targets*, **not proven universal acceptance limits**:
- Q01 controlled low-noise independent Airy: n_RMSE < 0.06, k_RMSE < 0.025, unseen-thickness R/T RMSE < 0.012 (automated).
- Synthetic low-noise cross-solver: forward-model parity should ideally be <0.005 R/T absolute RMSE **before inversion**; deviations require diagnosis, not relaxed thresholds.
- Target held-out RMSE must be selected from the measurement noise and independently justified by the forward-model discrepancy; use a fixed threshold declared before fitting.
- No numerical output counts as an identified n,k if multiple distinct branches are near optimal, if Jacobian effective rank is deficient, if parameters are pinned at bounds, or if unobserved-sample predictions fail. A single detected branch is *not* proof of uniqueness.

## 5. Competing branches protocol
Run at least 20 deterministic and randomized starts for research-grade ambiguity study across physically defensible n,k/dispersion/thickness parameterizations, keeping all solver options and observations identical. Preserve each start's initial state, termination code, evaluation count, final weighted objective, n/k probes and thickness. Cluster candidate curves on a declared wavelength probe grid; report pairwise n/k RMS distance and the number within 5% cost of the minimum (also examine noise-equivalent chi-square differences). Nonconverged starts are not evidence of uniqueness. Reference truth must never influence which branch is selected; evaluate accuracy only after ranking all candidates by the objective.

**Current implementation limitation:** multilayer fit runs 4 deterministic starts for splines and only 1–2 for other models; current branches use n/k RMS separation 0.05 and 5% relative cost. This is a preliminary diagnostic, NOT a comprehensive ambiguity search.

## 6. Held-out measurement protocol
Hold out entire independent specimen configurations, not randomly sampled adjacent wavelengths from the same trace. At least one of thickness, incidence angle, TE/TM polarization, layer order/known overlayer or substrate must change while shared-material n,k remain valid. The withheld configuration is not used to tune hyperparameters, initialization, model class, branch selection or bounds. Export predicted and measured R/T, residuals and metric tables. A wavelength-only holdout is explicitly weaker evidence due to spectral correlation.

## 7. Robustness and sensitivity
Run noise levels (0%, 0.2%, 1% absolute R/T standard deviation where physically meaningful), independent noise seeds, thickness perturbations ±5%, RI uncertainty of known layers, and reduced spectral bands. Summarize median and 95th percentile errors and solver failure rate. Do not automatically clip synthetic noisy data to [0,1] without recording the operation and its bias; prefer bounded measurement models or document the clipping.

## 8. Release gate
- All deterministic unit/API and browser checks pass on the exact deployed Git commit.
- Each core supported inversion family passes at least one independent analytic and one cross-solver or experimental reference test appropriate to its assumptions, with fixed acceptance criteria.
- Full report includes successes **and failures**, per-case figures and data, flags for competing branches and nonconvergence, provenance and reproducible run configurations.
- Never silently classify a fitted spectrum as accurately recovered n,k unless reference/identifiability diagnostics support it.
- Do not block a useful exploratory laboratory solely for ambiguity; label it 'exploratory / non-unique' and disable unsupported accuracy claims.

## 9. Current validated evidence
- Prior API, model-parity and energy/passivity regression tests pass on GitHub CI.
- `tests/test_multilayer_dispersion.py::test_independent_airy_blind_recovery_and_unseen_thickness` passed on GitHub Actions commit `ad6f5d97` with the Q01 thresholds above. Exact observed margins have not yet been exported into a machine-readable benchmark report.
- Q02–Q10, comprehensive multistart search, FMM cross-solver recovery and experimental validation **remain open**. Passing prior self-generated synthetic datasets must not be treated as completing them.
