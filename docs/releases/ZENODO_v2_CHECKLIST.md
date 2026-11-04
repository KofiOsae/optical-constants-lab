# Archival and DOI checklist for Optical Constants Laboratory v2.0.0

## Immutable target

The published GitHub release `v2.0.0` is tagged at exact commit:
`17c94963edd364f87e1ddddbe69c45d3c5c833d1`.

Do not move or recreate this tag. Current `main` contained that commit when this checklist was written. The `CITATION.cff` in the **archived v2.0.0 tag** still says version `1.0.0`; this is a known metadata defect. Future branches may correct the **current** citation file, but doing so does not modify the archived v2.0.0 bytes. In Zenodo's draft metadata, select **Software**, give release version **2.0.0**, and confirm the correct creators and custom license before publishing. Do not silently claim the v2 archive's embedded citation file was fixed.

## Required author/license verification before public DOI registration

- Verify creator names, order, affiliations, ORCID IDs (optional), and copyright ownership directly with the maintainer; do **not** infer the legal author from a GitHub account handle.
- The source license is **Research and Independent Non-Commercial Use License**, a custom source-available license, **not** OSI open source, MIT, GPL or CC-BY. Choose an appropriate `Other`/custom-rights option or record the exact license terms in Zenodo; do not select a falsely equivalent standard license.
- Confirm included third-party data, code, and images may be distributed in a permanent public archive. Exclude files lacking clear redistribution rights before publishing. Synthetic spectra from this repo are explicitly labeled as generated examples, not measurements.
- Check Zenodo draft preview for version, license, authors, keywords, checksum, public files, and the target code. Publishing creates a persistent public record.

## Preferred Zenodo path

1. Sign in to https://zenodo.org, connect GitHub in account settings, select GitHub from profile, sync, and enable `KofiOsae/optical-constants-lab`.
2. Inspect whether `v2.0.0` (published **before** enabling integration) can be archived from the Zenodo GitHub integration. Do not assume a historic release is ingested retroactively.
3. If it is not available, create a Zenodo **new upload** of the source archive downloaded directly from the existing GitHub `v2.0.0` tag. Verify file list and checksum before publishing; include a related identifier for the GitHub release URL and exact commit in the description.
4. Use title `Optical Constants Laboratory v2.0.0`, resource type `Software`, publication date as appropriate, release version `2.0.0`, and verified creator(s).
5. Use DOI automatically assigned by Zenodo upon publication. Record the **version DOI** and (if provided) **concept DOI** separately. Do not fabricate or preannounce an unissued DOI.
6. After DOI issuance, add the real DOI to `CITATION.cff`, README and application citation UI **in a new release**, without altering the existing `v2.0.0` tag. Never introduce a link to a draft DOI as if already registered.

## Suggested draft description

> Optical Constants Laboratory v2.0.0 is a browser-accessible research and teaching platform for retrieving wavelength-dependent optical properties from reflectance/transmittance measurements. It includes Single, Joint and Multilayer workflows, parametric and flexible dispersion representations, progressive and multistart fitting, diagnostic plots, reference comparisons, and data exports. The included teaching spectra are synthetic; numerical regression and cross-checks do not establish a unique inverse solution for arbitrary samples. The usage dashboard stores anonymous UTC-day event totals in a separate database. Repository version: v2.0.0. Exact release commit: 17c94963edd364f87e1ddddbe69c45d3c5c833d1.

Keywords: optical constants; thin films; reflectance; transmittance; inverse optics; dispersion; Kramers-Kronig; spectroscopy; scientific software.

## After publication

Capture the full Zenodo record URL/DOI, record whether it archives the GitHub tag or manually supplied source, and verify its code matches the exact release commit. Update the repository's archival log and citation references in a new development PR. Scientific data DOI(s) can be separate from a **software** DOI if future independently validated measurements are deposited.
