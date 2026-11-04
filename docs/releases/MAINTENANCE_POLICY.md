# Maintenance and release governance

**Stable production:** `main` → Render service `optical-constants-lab` (original public URL).  
**Test deployment:** `develop/v2-visual-publication` → Render `optical-constants-lab-v2-staging`.  
**Archived v2.0.0:** immutable tag `v2.0.0` → `17c94963edd364f87e1ddddbe69c45d3c5c833d1`.

## Change workflow

1. Branch new scientific or UI changes from `develop/v2-visual-publication`. Keep `main` and published tags untouched while developing.
2. Open PR to development. Add tests for model physics, API response shape, saved-project compatibility, browser layout and export files as applicable.
3. Require passing numerical regression, reconstruction and browser CI, plus new independent reference check. A green numerical check alone is not an assertion of experimental accuracy.
4. Deploy to staging and document manual acceptance with fixed test inputs for all three workflows, relevant model families, progressive fitting, multistart, charts, export and Save/Load. Check anonymous counters with Neon and verify no personal identifiers stored.
5. For release, compare development to production. Merge a reviewed release PR into `main` only when all gates pass. Render production auto-deploy follows main; verify the exact live commit, HTTP health, workflows and database event persistence.
6. Create a version tag on the **exact verified production commit**. Confirm tag SHA and release metadata. Never move version tags. Update citation metadata for *future* releases in the source before tagging.
7. Roll back only via a reviewed reversible commit or an explicitly selected earlier Render deployment. Preserve database schema compatibility; never delete production data casually.

## Mandatory release checklist

- [ ] version, release commit, date and changelog recorded
- [ ] PR review and all relevant workflows green
- [ ] synthetic forward reference verified independently of production implementation
- [ ] inverse identifiability limits/held-out and noise tests documented
- [ ] Single/Joint/Multilayer browser acceptance completed
- [ ] CSV/Excel/project and methods exports tested for supported models
- [ ] no secrets, personal visitor identifiers or uploaded spectra in usage database
- [ ] `CITATION.cff`, README, UI version and legal copyright/license reviewed
- [ ] production Render `live` commit equals intended release SHA
- [ ] DOI only published or cited once its live Zenodo record is verified

## Source and licensing note

The repo's license is a custom non-commercial source-available license, not OSI-approved. Do not describe it as unrestricted open source. Data and third-party licensing require separate confirmation for any permanent public archive.

## Release version policy

Use patch versions for backward-compatible bug fixes, minor versions for new backward-compatible capabilities, and major versions for breaking changes. Update associated scientific benchmark tables and validation results whenever an optical solver/model changes.

This document is operational guidance, not a claim that every historical release met all future gates.
