# Optical Constants Laboratory v2 — browser acceptance and research validation
Version: staging, **not** a v2.0.0 release. Test each feature in staging and attach screenshots or exported project files. Never overwrite the published v1.0.0 tag.

## Baseline and user interface
- [ ] Three workflows visible, live stack diagrams for Single/Joint/Multilayer, diagrams collapse completely, controls focus corresponding layer.
- [ ] Main Plotly plots remain interactive and usable (zoom, hover, legend).
- [ ] One optional Save plot control and one Publication & reproducibility control per active view; no duplicate controls or large empty page regions.
- [ ] Export PNG and SVG plot; numerical plotted-data CSV has correct headers/units and curve values.
- [ ] Export methods TXT, Markdown, and browser Save as PDF, and portable reproducibility JSON. Reopen original project using the standard project-import control (the reproducibility JSON is a wrapper; use its `project` field).
- [ ] Verify that every export includes assumptions and does not invent a DOI, solver commit, or verification status.

## Physics exercise A: transparent single film
File: `transparent_film_120nm.csv`. Input n=1.65, k=0, d=120 nm, ns=1.46, normal incidence. Use transparent model, thickness fixed, mixing 0, noise 0. Check R+T ~ 1. Compare R/T fit and recovered n,k to reference, try several initializations and document ambiguity.

## Physics exercise B: transparent joint films
Files `transparent_film_120nm.csv`, `transparent_film_240nm.csv`. Thicknesses fixed at 120 and 240 nm respectively. Shared unknown n=1.65,k=0. Compare per-film reconstructions and common fitted n,k. Remove a spectrum and test sensitivity.

## Physics exercise C: absorbing film
File `absorbing_film_180nm.csv`. d=180 nm, ns=1.46; n(λ)=1.73+0.12(550/λ−1), k(λ)=0.07 exp[-(λ−610)^2/(2·60^2)] nm-convention. Fit cautiously with exploratory flexible model. This is not a globally causal material reference; do not treat failure of KK-based fitting as a solver failure.

## Physics exercise D: multilayer with cap
File `multilayer_cap50_film150.csv`: Air / **known 50 nm, n=1.8 cap** / **unknown 150 nm, n=1.65, k=0 film** / n=1.46 substrate. Verify order in the stack diagram. Check the effect of accidentally switching the known layer to `after` then restore correct order. Fit fixed unknown thickness, zero mixing/noise.

## Automated reference tests
`tests/test_v2_examples.py` compares exact synthetic R/T files to the laboratory's coherent-stack solver, checks passivity and transparent energy conservation. Run with `python -m pytest -q tests/test_v2_examples.py`. This verifies forward-model consistency, **not** parameter-recovery identifiability or independence of the reference.

## Required evidence per experiment
Record: environment and deployed commit, source CSV, layer order, substrate assumptions, selected model, optimizer settings/bounds, fitted parameters, plot exports, R/T residual metrics, ground-truth n,k errors, and whether a second optimizer initialization agrees. Label unsupported and non-identifiable cases explicitly.

## Release gate
Do not merge to production or tag v2.0.0 before: (1) CI passes; (2) browser tests pass in all three workflows; (3) each example is reproducible; (4) all fit/forward-model limitations are documented; (5) published figures and TXT/PDF/JSON contain honest version and provenance metadata.
