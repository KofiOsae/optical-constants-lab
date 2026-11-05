# Optical Constants Laboratory — Browser / Tablet / Render v0.3

## Official scientific software citation (v2.0.0)

Asamoah, Benjamin O. (2026). *Optical Constants Laboratory* (Version 2.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23966359

- [Version-specific Zenodo archive](https://doi.org/10.5281/zenodo.23966359)
- [Immutable GitHub v2.0.0 release](https://github.com/KofiOsae/optical-constants-lab/releases/tag/v2.0.0) — commit `17c94963edd364f87e1ddddbe69c45d3c5c833d1`
- [Public laboratory](https://optical-constants-lab.onrender.com)

**Archive note:** The original v2.0.0 GitHub source archive is immutable and predates DOI registration; this README citation update is a later development commit and does not modify the version-tagged release.

Same scientific Python engine with responsive HTML frontend. Run locally on a computer or deploy the Flask API on Render using `render.yaml`. A tablet accesses the hosted site through its ordinary browser; it does not run Python itself.

**Quick start (computer):** install `requirements.txt`, then `python server.py` and open `http://127.0.0.1:8765`. **Tablet:** follow `DEPLOY.md` to deploy and open the hosted service URL.

This is a deployment-ready *configuration*, not a confirmed cloud deployment. The frontend, engine, and four synthetic datasets are retained from v0.2. See `TUTORIAL.md` and `DEPLOY.md` for scientific scope and deployment cautions.

### Custom experimental spectra (v0.3.7)
Upload CSV/TSV/TXT/JSON through **Custom experimental spectra**. Use columns `wavelength_nm` or `wavelength_um` and `R` and/or `T`. Choose whether signals are fractions or percent; enter known substrate index and film thickness estimate. The server runs the existing spline inversion without accessing R6G ground truth. This first importer supports a single coherent film on a semi-infinite nonabsorbing substrate. An initial film thickness is mandatory and can optionally be fitted within ±15%. The custom upload API is `/api/custom-fit`. FMM data export may be converted to this common tabular schema; automatic FMM JSON import is not yet supported.

## v0.3.8: General spectral import preview

The Custom experimental spectra panel now previews uploaded wavelength-indexed n/k, R/T, ellipsometry Ψ/Δ, absorbance, absorptance, and dielectric-function tables without requiring the optional R6G example. The `/api/import-preview` endpoint validates files locally on the app server and returns detected columns, converted wavelengths, and data for plotting. It does not perform inversion for n/k or ellipsometry; the existing R/T inverse solver remains unchanged. Ambiguous wavelength units must be explicitly selected. R6G can be loaded on demand via a separate button, and its failure does not block uploaded data.


## Copyright, research use and commercial permission

Copyright © 2026 Optical Constants Laboratory contributors. The original Laboratory source code and documentation are made available at no license fee for personal, educational, academic and independent **non-commercial** research, under the conditions in [LICENSE](LICENSE). Commercial use (including paid services and proprietary integrations) requires separate written permission from the relevant copyright holder(s), requested through the repository maintainers. Third-party dependencies and material data retain their own licensing terms. This is a custom source-available non-commercial license, not an OSI-approved open-source license.
