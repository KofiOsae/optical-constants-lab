# Independent numerical validation protocol

Version v1.0.0 is released software, not a claim of independently verified physical accuracy. Independent researchers are invited to run these tests and publish code, data, deviations, and negative findings.

## Reproduction
1. Check out immutable tag `v1.0.0` (commit `57441470f272025b020e36bfd4d29b24f897ba73`).
2. Record Python/library versions, OS, spectral units, incidence angle, polarization, layer order, noise, optimizer settings, initializations, bounds, wavelength coverage, and numerical tolerances.
3. Run the repository's CI and analytical benchmarks and report any failures separately from browser UI issues.

## Independent experiments
- **Transparent single film:** Generate Fresnel/TMM R and T with independent code, recover a known dispersion and thickness; vary interference-fringe density and noise.
- **Joint films:** Use two or more distinct known thicknesses of one material, fit the common n,k and report both reconstruction RMSE and parameter error.
- **Absorbing substrate:** Independently evaluate R/T for a film over complex `n+ik`, including energy-flux conventions.
- **Coherent finite substrate:** Independently compute front- and back-interface interference and compare spectral phase-sensitive response. Do not treat this as validation for incoherent thick plates.
- **Effective interface mixing:** Create a known optical interlayer, then compare fitted effective thickness and n,k bias against a zero-mixing control. This is not an RMS roughness or diffuse-scattering validation.
- **KK/SSKK:** Check analytic causal Lorentz dielectric functions with independently generated broad-band references; quantify finite-band truncation, sensitivity to the anchor and epsilon-infinity, and interior/full-band residuals.
- **Model identifiability:** Test distinct initializations, oscillator orders, thickness mismatch, noise and held-out spectra. Publish Jacobian rank, condition number, and alternate minima.

## Reporting
Supply the raw R/T spectra, input optical constants, generation solver, wavelength grid, fitted n,k, residuals, optimized parameters, exact commit, environment and executable scripts. For each test, state predeclared numerical acceptance tolerances justified by precision, wavelength sampling and noise, then report pass/fail. Include discrepancies rather than selecting only successful fits.

## Release caveat
Initial user acceptance tests covered saved projects, batch assignment, interactive layout, and KK/SSKK display. The final independent reconstruction-regression acceptance test remained unverified when v1.0.0 was published. This protocol supports subsequent external verification without changing the archived release.
