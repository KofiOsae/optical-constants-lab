# Original R6G Workbench dataset: verification record

This package embeds the **actual** archived Workbench spectrum (`data/original_plot_raw_7.json`, 201 wavelengths from 300 to 1000 nm) and reference material (`data/R6G_Singh.csv`). It also embeds the previous Phase I fit archive (`data/phase1_fit_verified.npz`).

## Verified here

- Workbench model: air / uniform R6G film (1.000 µm) / air, normal incidence.
- Original exported R, T and A sum to unity to machine precision.
- Independent Python coherent Fresnel forward model using tabulated R6G n/k reproduces the exported R and T within the numerical regression threshold of 0.005 RMSE each.
- Previous Phase I fit archive records a thickness of 999.37988 nm and contains the same R/T input samples as the Workbench export.
- All 12 automated tests pass in this environment, including four new original-dataset regression checks.

## Not yet verified

- The new browser application's `/api/fit` endpoint has **not** reproduced the original Phase I 999.38 nm inversion. Its current spline fit is a different inverse parameterization and lacks the original Phase I algorithm. The archive's 999.38 nm result is a historical benchmark, **not a new inversion result**.
- The web server and tablet UI were not tested end-to-end (Flask is unavailable in this execution environment).
- No new Render service has been deployed from this package.

Run `PYTHONPATH=. pytest -q` in the project root after installing requirements.

## v0.3.6 reference-loading correction
The bundled R6G CSV header is `wavelength_um,n,k`, but the prior `/api/reference` parser only recognized `wavelength`, `wl`, or `lambda`. Consequently the reference endpoint raised `Missing column` and the browser never loaded the examples. The parser now explicitly recognizes `wavelength_um` and `wavelength_nm`, converts µm to nm, and checks finite and monotonic data. Added regression tests for the original file and loader. Full HTTP/browser verification and online deployment remain outstanding because Flask is not installed in this execution environment and the standalone service has not been deployed.

## v0.3.9 upload pipeline hotfix

- Custom R/T fitting now calls the same `inspect_material` validator used by import preview. Both use identical wavelength and signal units, column aliases, sorting, and validation.
- `R_percent` and `T_percent` are accepted as percentage-valued input columns.
- Fit button stays disabled until a compatible R/T preview succeeds; a new file selection invalidates the previous preview.
- API calls surface network failures and non-JSON HTTP errors rather than generic JSON parse errors.
- Flask returns JSON for 413 upload-too-large responses.
- Tests: `python -m pytest -q tests/test_material_import.py tests/test_shared_upload_pipeline.py tests/test_r6g_reference_units.py` (9 passed in packaging environment). Full suite requires Flask installed.
- Browser JS syntax checked with Node.js. Render deployment not yet performed.
