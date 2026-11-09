# Nanoparticle Optics → Optical Constants Laboratory integration contract

Status: approved design requirements; **not** a claim of completed Meep/Colab implementation.

## Final destination
Integrate as a dedicated **Nanoparticle Optics** tab in Optical Constants Laboratory after research-only engineering and scientific validation. Preserve existing Single, Joint, Multilayer, and published v2.0.0 functionality. The research preview remains isolated until approval.

## User flow
1. Configure material n,k or epsilon1,epsilon2 spectra (with provenance), host, particle geometry, 2D/3D, boundary conditions, source, resolution, monitors, wavelength range, simulation time, and expected outputs.
2. Validate geometry and physics. Explicitly distinguish 2D infinite-cylinder from finite 3D particles; estimate memory and computational cost for 3D, show warnings and numerical-convergence recommendations.
3. Create a **self-contained Colab notebook** with all geometry, dispersive materials, simulation inputs, and visualization/export code embedded. Notebook must not require manual reconstruction of inputs. Version the job schema and notebook generator.
4. Provide a direct **Open in Colab** action for a notebook reachable through an authorized public or user-controlled storage URL. Do not claim that a local notebook can be launched through Colab without an accessible URL. A one-click notebook download plus Colab upload is a fallback, not the target.
5. User connects a Colab runtime and clicks **Run all**. The notebook installs or verifies a compatible Meep runtime, runs the model, plots fields and spectra, and exports CSV/NPZ with machine-readable metadata. Colab execution is user-initiated; no silent remote execution or automatic access to private Drive.
6. Provide **Import Meep Results** in the Laboratory. Validate schema, units, simulation dimension, normalization and source/monitor metadata. Plot and compare against Mie, TMM, ellipsometry and experimental spectra. Permit subsequent optical-constants inversion only where the measured/simulated observable and model justify it.

## Solver and science constraints
- 2D: invariant cylinder and other explicitly supported planar geometries. 3D: finite sphere/ellipsoid or explicitly meshed geometry; never equate 2D per-length scattering with 3D cross section.
- Dispersive complex epsilon requires Meep-supported passive material models (e.g., fitted Lorentz/Drude susceptibilities) or an explicitly limited narrowband approximation. A CSV spectrum cannot simply be supplied as a frequency-independent complex epsilon for broadband FDTD.
- For scattered-field cross sections, calculate incident normalization and scattered flux with reference/background subtraction. Document 2D units versus 3D units.
- Track PML, resolution, time decay, convergence and energy balance; flag nonconverged results.
- Results must include schema version, job ID, dimension, geometry, wavelength_nm, polarization, units, incident normalization, solver version, material provenance, and monitor definitions. Optional arrays: R,T,A, scattering, absorption, field samples.
- Colab free-tier resources and Meep package availability are variable; 3D runs may require smaller jobs or a paid/local runtime. Validate installation in the generated notebook and report failures rather than promise guaranteed execution.
- No unapproved publishing of uploaded proprietary spectra or data; embed user inputs in the notebook only with clear export/privacy disclosure.

## Acceptance criteria
- Notebook generator has deterministic tests for both 2D and 3D.
- Colab notebook opens with prepopulated geometry and materials and requires no manual re-entry; user only authorizes the runtime and clicks Run all.
- 2D and 3D smoke tests complete on a supported runtime and export round-trippable results.
- Import rejects unit mismatches, missing normalization, and incompatible 2D/3D comparisons.
- Validate against analytic cylinder/sphere cases and convergence sweeps, followed by literature-backed scientific validation.
