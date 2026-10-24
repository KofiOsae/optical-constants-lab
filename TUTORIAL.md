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
