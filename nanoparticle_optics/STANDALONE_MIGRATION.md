# Nanoparticle Optics Workbench — standalone migration

This research application is deployed independently from Optical Constants Laboratory v2.0.0.

## Current isolation
- Flask entry point: `nanoparticle_optics.preview_server:app`
- Python package: `nanoparticle_optics/`
- Dedicated Render service: `nanoparticle-optics-workbench`
- Research branch: `research/nanoparticle-optics-ema`
- Automatic deployments disabled until engineering acceptance.
- No production Optical Constants Laboratory routes or releases are modified.

## Repository extraction (pending)
Create a dedicated `nanoparticle-optics-workbench` repository. Copy the
`nanoparticle_optics/` package, `tests/test_nanoparticle_*.py`, workflow,
and licensing/provenance documentation. Adjust package paths and deploy config,
then switch the dedicated Render service to the new repository after tests pass.

## Cross-workbench interoperability contract v0.1
Primary material CSV:
```csv
wavelength_nm,n,k
500,1.50,0.00
600,1.49,0.00
```
Permittivity alternative:
```csv
wavelength_nm,epsilon1,epsilon2
500,2.25,0
600,2.2201,0
```
Conventions: vacuum wavelength in nm; passive exp(-i omega t) convention,
`n+ik` with `k>=0`, `epsilon=epsilon1+i epsilon2`.
No wavelength extrapolation. Specify whether values are intrinsic particle,
host, or effective-film optical constants. Preserve provenance, temperature,
measurement method, uncertainty, and applicable spectral range in an adjacent
metadata JSON file where available.

Output spectra: `wavelength_nm,R,T,A` are fractions, not percentages.
Document coherent/specular versus hemispherical/integrating-sphere measurements.
Do not silently interpret coherent film absorption as diffuse scattering.

FMM interchange: material spectra only until geometry and polarization
conventions are explicitly agreed. Meep interchange: geometry-specific field
slices with spatial units, frequency, resolution, PML, source and solver
metadata; do not interpret field snapshots as normalized cross sections.

## Engineering acceptance
1. Execute Python regression tests and route checks.
2. Verify all navigation modes and plot observables, including axes and legends.
3. Verify material imports, fill-factor comparison and pore infiltration.
4. Verify Meep notebook on a fresh authorized Colab runtime.
5. Confirm downloadable files and metadata are reproducible.
6. Keep scientific/literature validation and ellipsometry stage deferred.
