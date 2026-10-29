# Optical Constants Laboratory — Theory and Guided Tests

## 1. Scope and status
This prototype provides a coherent single-film Fresnel forward model and a regularized B-spline inverse model for four synthetic R6G cases. The supplied R6G optical constants are a reference *for numerical recovery*, not independent experimental truth. The following methods are **not implemented** in this version: original Swanepoel envelope extraction, causal Lorentz/KK inversion, measurement-file upload, calibrated uncertainty intervals, and automated measurement design.

## 2. Definitions and sign convention
Use the time convention exp(-i omega t) and complex index N=n+ik with k>=0 for passive materials. The absorption coefficient is alpha=4 pi k/lambda. Single-pass intensity attenuation is exp(-alpha d), which is not the full thin-film transmittance because it neglects interfaces and interference.

## 3. Coherent thin-film model
For ambient index n0, film complex index N, and semi-infinite substrate index ns, define r01=(n0-N)/(n0+N), r12=(N-ns)/(N+ns), delta=2 pi N d/lambda. At normal incidence,

    r=(r01+r12 exp(2i delta))/(1+r01 r12 exp(2i delta))
    t=(t01 t12 exp(i delta))/(1+r01 r12 exp(2i delta))
    R=|r|^2,  T=Re(ns)/Re(n0) |t|^2

for nonabsorbing entrance and exit media. The application also models oblique incidence and s/p polarization. A finite glass slab with backside reflection needs a different (often mixed coherent–incoherent) model; the built-in glass substrate is semi-infinite.

## 4. Original Swanepoel method (theory only)
For a suitable film on transparent substrate, upper and lower transmission envelopes TM and Tm yield an approximate film index via

    A=2 s (TM-Tm)/(TM Tm)+(s^2+1)/2
    n=sqrt(A+sqrt(A^2-s^2)).

Thickness follows from interference order/adjacent extrema, subject to consistent wavelength ordering. This approach requires resolved fringes, appropriate substrate assumptions, and weak/moderate absorption. **It is not an executable extraction option in this prototype.**

## 5. Why B-splines?
The forward model remains Fresnel theory. A spline only parameterizes the unknown spectra:

    n(lambda)=sum_j a_j B_j(lambda)
    k(lambda)=sum_j b_j B_j(lambda).

Instead of fitting independent n,k at every wavelength, one optimizes a limited set of coefficients. Smoothness can be encouraged by penalizing curvature. A smooth fit does **not** establish Kramers–Kronig causality or unique parameter recovery. More basis functions can reduce fit error while increasing instability.

## 6. Multiple-spectrum fitting
The four synthetic cases are air/R6G(1000 nm)/air, air/R6G(50 nm)/air, air/R6G(1000 nm)/glass, and air/R6G(50 nm)/glass. The same n(lambda),k(lambda) describe all cases, while sample thickness and substrate differ. Joint fitting uses a common material model, not independent fitted indices averaged afterward. It is only justified experimentally if the samples have the same intrinsic material response.

## 7. Optimization methods
Least squares is a local solver; multistart repeats local fits; differential evolution and Harris Hawks Optimization (HHO) explore parameter space globally; HHO-hybrid follows global search with local refinement. None guarantees physically correct optical constants when the data are non-identifying. Compare methods with matched forward-evaluation budgets, multiple random seeds, and withheld measurements. HHO is an *optional addition*, not the default.

## 8. How to interpret results
Inspect (i) reflectance and transmittance residuals, (ii) optical-constant error against known synthetic reference, (iii) thickness recovery, (iv) optimizer success and evaluation count, and (v) stability under basis complexity and different starting guesses. The lowest spectral RMSE alone does not prove a unique or causal reconstruction. Error bars derived from a local Jacobian can be misleading for multimodal fits.

## 9. Guided experiments
**A. Forward physics:** Compare thick and thin R6G near its absorption band. Calculate exp(-4 pi k d/lambda) to anticipate transmission suppression; compare with the complete Fresnel result.

**B. Substrate effect:** Compare the same thickness with air and glass exit media. Explain differences in interface reflectance and interference, without interpreting them as material-index changes.

**C. Fit one versus four spectra:** Use least squares, a fixed basis size, and the same iteration budget. Compare RMSE of n and k, not only spectral RMSE. Repeat with unknown thickness if supported by the interface.

**D. Solver comparison:** Compare least squares, multistart, differential evolution, HHO, and hybrid HHO using identical data and bounds. Stochastic results require repeated seeds and matched computational budgets before declaring a winner.

**E. Complexity:** Increase spline basis count and adjust regularization. Observe whether spectral residual decreases while recovered optical constants become less stable.

## 10. Verification versus validation
**Verification** checks equations and code: zero-thickness and index-matched limits, R+T=1 for lossless films, convergence and independent solver agreement. **Synthetic validation** checks recovery against known n,k,d under noise. **Experimental validation** additionally requires calibrated measured R/T, independent film thickness, measured substrate constants, and an accounting of diffuse scattering and backside reflections. Passing synthetic tests does not establish experimental accuracy.


## 11. Multilayer laboratory: selecting an inversion model (October 2026)
The current multilayer workspace fits **one shared unknown isotropic material** across 1–12 measured R/T spectra with different known layers and optional thickness bounds. It uses coherent finite layers and a semi-infinite lossless substrate; do not interpret a fitted result as valid for incoherent substrate-backside reflections, diffuse scattering, graded layers, or anisotropy.

### Available inverse models
| Method | Use when | Assumptions and pitfalls |
|---|---|---|
| B-spline | Material identity or dispersion shape is unknown; exploratory reconstruction | Flexible but **not automatically causal**. Vary basis count and smoothness; check competing branches and sensitivity. |
| Lorentz | Bound-electron resonances dominate | Passive oscillators, but resonance number and initial guesses matter. |
| Drude–Lorentz | Free carriers plus interband/bound resonances may be present | Test multiple oscillator counts; a lower residual alone does not validate the recovered constants. |
| KK-constrained fixed-pole Lorentz basis | You want a causal, passive basis with fixed pole locations | This implementation is a restricted oscillator basis, **not** a general model-free KK inversion; inadequate poles can bias recovery. |

### Diagnostic-only options
Finite-band **KK** and singly subtractive **SSKK** checks assess a reconstructed spectrum; they **do not constrain the optimizer**. Finite bandwidth can produce truncation errors. SSKK requires an independently justified anchor.

### Other models and methods (not currently executable in the multilayer inverse)
- **Cauchy:** transparent, weakly dispersive dielectrics over a limited wavelength band; not suitable across absorption features.
- **Sellmeier:** transparent dielectrics away from resonances; do not extend a transparent-band fit through absorption.
- **Swanepoel envelopes:** preliminary index/thickness estimates for transparent or weakly absorbing films with resolved transmission fringes and suitable substrate conditions.
- **Drude-only:** free-carrier response when bound-electron terms can be neglected.
- **Tauc–Lorentz:** some amorphous semiconductors with a meaningful optical band edge; needs appropriate causal implementation and bandgap constraints.

### Practical decision workflow
1. **Define the optical stack:** illumination side, order of layers, thickness estimates, known-layer n,k, substrate index, incidence angle, polarization and measurement units. An unknown layer is not automatically equivalent to a known reference material.
2. **Inspect measured R/T:** check 0–1 fractions, spectral overlap, calibration, R+T, interference visibility, absorption and uncertainty. R+T<1 can arise from absorption or unmodeled scattering, so do not assume every loss is intrinsic k.
3. **Start reference-independent:** fit a low-complexity B-spline when the material is unknown, and test Lorentz/Drude–Lorentz where physically plausible. Do **not** seed with a reference n,k spectrum in a blind-recovery experiment.
4. **Use physically justified bounds:** n>0 and k>=0 for passive media, plus independent thickness limits. If parameters hit a bound, widen it only when supported by the physics; report bound sensitivity.
5. **Compare candidates:** inspect R/T residual structure, fit convergence, distinct near-optimal optical-constant branches, Jacobian rank/conditioning, and model complexity. A low RMSE or one converged start is **not proof of uniqueness**.
6. **Validate predictive power:** reserve independent spectra or configurations for held-out prediction when possible; perturb thickness, model complexity and starting points. A reference n,k file, if supplied, is used **only after fitting** for optional comparison, not for fitting.
7. **Report limitations:** distinguish forward-model agreement, inverse reconstruction, identifiability and experimental validity. If branches compete, report ambiguity rather than a single definitive n,k.

### ML-01 benchmark interpretation
For air / Au (30 nm) / silica (n=1.46), 400–1000 nm, a separately uploaded Au reference produced forward-only combined R/T RMSE ~0.001972 against FMM-generated R/T. This verifies approximate forward-model agreement for that configuration, **not** a successful blind inversion. Continue ML-01 using measured R/T and independently specified geometry only. The reference is an evaluation benchmark and must not enter optimization.


## Organic dyes: expanded causal oscillator models (development branch)
The planned additional models use a common causal complex dielectric response and the same
R/T forward geometry as Lorentz fitting. The existing **model dropdown** is the intended
point of access; no additional reconstruction workspace is required.

* **Modified Lorentz** uses a nonnegative combination of narrow and broad Lorentz
  responses at one center energy. Set strength, center, base damping, and broad fraction.
* **KK-consistent Gaussian broadening** averages passive damped oscillators over a
  distribution of resonance energies. Set strength, center and Gaussian sigma;
  homogeneous damping is small and prescribed in this experimental definition.
* **Voigt-type broadening** averages homogeneous Lorentz poles over a Gaussian
  energy distribution. Set strength, center, Lorentz gamma and Gaussian sigma.
* **Vibronic multiplet** uses a finite Huang–Rhys distribution of molecular transition
  sidebands. Set strength, origin energy, vibronic spacing, homogeneous damping, and S.

The dielectric convention is N = sqrt(epsilon) = n + i k with exp(-i omega t), and
positive oscillator strengths, energies and widths. Do not fit an arbitrary Gaussian
to k independently of n. Avoid assuming the number of fitted peaks establishes
a unique microscopic transition assignment.

**Worked synthetic protocol.** Construct a 450–950 nm reference spectrum for a 150 nm
film on a known n=1.5 substrate using one Voigt-like oscillator with
epsilon_inf=2.5, strength=0.8, center=2.3 eV, gamma=0.12 eV and
sigma=0.15 eV. Simulate R and T through the existing coherent forward solver.
Fit the generated data at fixed thickness and compare n, k, R, T and parameter
recovery, including 470–600 nm residuals; then perturb with a recorded noise seed
and repeat for a second thickness. Compare with 1–4 Lorentz oscillators and
held-out spectra. Values are **illustrative**, not published dye properties.

**Release caution:** This tutorial paragraph documents the science and intended
controls; it must not be interpreted as proof that the current development
branch has full UI/model-export coverage. Consult the release checklist.


## Interpreting fit diagnostics and oscillator initialization

The fit's **weighted RMSE** measures the normalized differences between modeled and measured reflectance/transmittance. It is not the same as an absolute R/T fraction or the optical-constant error. Inspect both per-film R/T RMSE and an independently uploaded reference n,k when available. A successful optimizer termination means a numerical stopping criterion was satisfied; it does **not** establish physical accuracy.

For a residual vector r and fitted parameter vector p, the local Jacobian is J_ij = partial r_i / partial p_j. Singular-value decomposition yields singular values sigma_1 >= ... >= sigma_m. The reported **condition number** is kappa(J) = sigma_max / sigma_min, provided the smallest singular value is positive. A full-rank Jacobian means parameters are locally distinguishable in the linearized fit, but neither full rank nor a modest condition number proves global uniqueness. Condition numbers depend on parameter units, scaling, chosen wavelengths, noise weighting and model. There is no universal numerical pass threshold. Weak directions may exhibit large parameter uncertainty even when the optical R/T fit appears good.

**Worked comparison from synthetic two-film tests:** Lorentz produced kappa ~11.8; Voigt ~17.2; modified Lorentz ~88.3; refined KK-Gaussian ~57.0. These values describe the local parameter scaling of the specific tests, not a ranking of material models. An earlier broad KK-Gaussian fit had kappa ~364 and conspicuous narrow numerical peaks, requiring separate quadrature-convergence checks.

**Choosing initial parameters:** epsilon_inf controls the long-wavelength dielectric baseline; oscillator strength controls amplitude; E0 sets the transition energy in eV (E = 1239.841984 / wavelength_nm); Lorentz gamma is homogeneous damping; Gaussian sigma is inhomogeneous energy-distribution width. The modified-Lorentz broad fraction mixes two linewidths, and the vibronic spacing and Huang–Rhys S control sideband positions and relative strengths. To demonstrate genuine recovery, rerun a known synthetic benchmark from deliberately perturbed initial guesses and physically meaningful bounds, not only from the known answer. A bound-clamped parameter, high correlation or strong dependence on initial guesses warrants a caution.

**Recognizing systematic shifts:** Compare fitted and reference n,k at matching wavelengths; inspect signed residuals, mean offsets and wavelength dependence. Incorrect substrate refractive index can produce a convincing R/T fit with biased n,k: in the controlled two-film Lorentz test, using substrate n=1.46 instead of the true 1.50 produced nonzero reference errors even though optimization converged. Verify film thickness, incidence angle, polarization, spectral units and input R/T channel calibration before changing dispersion physics. The application may flag offsets but cannot infer their unique cause from residuals alone.

**Interface scattering beyond the current mixing approximation — not yet available:** physical scattering roughness is not generally equivalent to a homogeneous effective-medium interfacial layer. Any future roughness control must specify whether it implements interfacial mixing, coherent reflection attenuation, or a scattering model, and its limits must be validated independently for single-film, joint-film and multilayer stacks. Keep roughness disabled for established baseline benchmarks.


## Effective interface mixing (roughness approximation)

The optional **interface mixing thickness** in nm models a thin, homogeneous, 50:50 mixture of the complex dielectric functions on each side of a boundary:

    epsilon_mix(lambda) = [epsilon_left(lambda) + epsilon_right(lambda)]/2
    N_mix(lambda) = sqrt(epsilon_mix(lambda))

The inserted finite mixing layer is propagated using the same coherent transfer-matrix optics as the rest of the stack. This is a **phenomenological interdiffusion / index-grading approximation**, not a direct measurement of RMS surface-height roughness or diffuse scattering. In particular, it does not predict roughness-scattered power, lateral correlations, haze, or polarization conversion. A thickness value describes the effective layer added at an interface; it should not be identified automatically with an AFM-measured Sa or Sq value.

**Single spectrum:** use the effective interface mixing field (nm); zero disables it. **Joint spectra:** a common interface mixing thickness is applied to each measured film. **Multilayer:** enter a comma-separated value for every boundary, in order from ambient through all finite layers to substrate. For a film with no additional known layers, enter two numbers, e.g. `2, 3`, for ambient/film and film/substrate; blank leaves all interfaces ideal. With one known layer above and one below, enter four values. The smooth-interface result must be unchanged when all values are zero.

Start with 0 nm, then compare fixed 1–5 nm hypotheses under the same measured spectra and material assumptions. Look for lower held-out R/T errors as well as stable n,k and thicknesses. Effective interface mixing can correlate with bulk absorption and film thickness: a lower training RMSE alone does not establish real roughness. Avoid making broad scattering claims without angular or diffuse measurements.


## Choosing a substrate and measurement geometry before inversion

**Applicability gate (current implementation):** The laboratory's coherent film and stack solvers require an effectively semi-infinite, nonabsorbing ambient and substrate in the modeled wavelength range; finite layers may absorb. The laboratory **does not currently model a thick substrate's rear interface or incoherent internal reflections**. In this context, 'semi-infinite' means the rear surface has negligible influence on the collected signal, not that real glass is literally infinitely thick. Ordinary transparent slides can have significant backside reflections and interference; verify the measurement configuration before using this approximation.

For measured samples, select a flat, optically polished, chemically compatible substrate with reliably known refractive index n_s(lambda), low absorption (k_s approximately 0) throughout the **entire measured band**, and negligible unmodeled substrate photoluminescence, scatter, or backside reflection. Candidate material families include optical glass, fused silica/quartz, and sapphire, **only within their respective verified transparency windows**. Avoid assuming n_s=1.50 for every piece of glass; use substrate-specific calibrated data where possible. Beware birefringence (e.g. sapphire), dispersion, fluorescence and coating-induced substrate changes. A blank-substrate R/T measurement is strongly recommended using the same optics, beam footprint, detector settings, polarization and incidence angle.

**Before uploading experimental R/T:** record wavelength range, substrate manufacturer/type, n_s(lambda), thickness, film deposition side, film thickness estimate and uncertainty, polarization, angle of incidence, reference standard and whether specular or total R/T was collected. Check whether the instrument separately removes the blank-substrate signal, corrects for multiple reflections or captures the rear surface; do not automatically treat normalized or substrate-referenced data as absolute sample R/T. A substrate correction is not equivalent to a physically matched forward model.

**Current limits:** weakly or strongly absorbing substrates require complex N_s=n_s+i k_s and appropriate reflectance/transmittance boundary conditions; the current kernels reject absorbing substrate indices. Opaque substrates require reflectance-only inversion with the corresponding complex substrate optics. Finite transparent or absorbing plates may need incoherent, coherent or hybrid treatment depending on thickness, spectral resolution and coherence. An unmodeled substrate can bias retrieved film n,k, thickness and inferred interface mixing.

**Recommended decision:** If the measured sample, blank-substrate response and backside configuration cannot be represented by the current model, do not treat the extracted n,k as quantitatively validated. Instead, use documented **synthetic R/T spectra generated with the laboratory's matching forward assumptions** for numerical solver validation. Matching synthetic data establishes algorithmic consistency, not experimental applicability. Report model suitability and uncertainty separately.

**Controlled substrate/roughness lesson (two Voigt films, 150 and 300 nm):** with true substrate n_s=1.50 and mixing 0 nm, the known synthetic optical constants were reproduced to numerical precision (weighted RMSE 6.716e-14). With correct n_s and 3 nm mixing, the weighted RMSE became 0.728456 and reference n RMSE 0.009416. With an incorrect substrate value, 3 nm mixing reduced R/T error while making n,k recovery worse: extra interface parameters can compensate for incorrect geometry. Do not equate a lower fitting loss with more physically accurate optical constants.

