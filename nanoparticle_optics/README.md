# Nanoparticle optics: isolated research extension

Branch: `research/nanoparticle-optics-ema`. No production routes, imports,
requirements, database migrations, or deployment configuration are modified.

## Scientific scope
The first implementation provides spherical-inclusion Maxwell-Garnett and
symmetric two-phase Bruggeman effective permittivities. The time convention
is exp(-i omega t), with passive Im(epsilon)>=0. The models are quasistatic:
particle diameter must be much smaller than wavelength *in the host*.
They do not describe size-dependent scattering, aggregation, percolation
physics, anisotropy, or near-field interactions.

The Bruggeman quadratic has two roots; the initial root selection uses
proximity to the Maxwell-Garnett estimate. This is **not** a fully validated
continuation algorithm near resonances or high-contrast transitions.
Before release, implement spectral branch tracking, passivity checks,
independent numerical comparisons, and failure diagnostics.

## Validation status
Tests supplied but NOT executed in the connected GitHub runtime. Status:
**unverified / research only**. No integration or deployment approval.

## Milestones and gates
1. Repository-wide release/architecture/license/test inventory; pin baseline.
2. EMA: analytic limits, complex-valued benchmarks, continuity/passivity.
3. Mie homogeneous spheres: cross-sections, Rayleigh limit, optical theorem,
   absorption nonnegativity, established independent implementation.
4. Size distributions: number/volume weighting, normalization, convergence.
5. Material metadata: schema, units, provenance, licenses, import/export.
6. Ellipsometry inverse: synthetic recovery, identifiability and uncertainty.
7. Optional external Meep 2D: explicit opt-in, separate execution, explain
   infinite-cylinder / 2D-geometry assumptions vs 3D spheres.
8. FMM interchange: documented schema and round-trip fixtures.
9. Full regression, scientific report, manual approval; only then integration.

No third-party code or data has been copied into this extension.
