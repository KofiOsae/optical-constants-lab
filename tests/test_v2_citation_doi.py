"""Regression guard for versioned Zenodo software citation references.

The published v2.0.0 Git tag is immutable; this validates development
metadata and in-app citation text, not the historical tag contents.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOI = "10.5281/zenodo.23966359"
RELEASE_SHA = "17c94963edd364f87e1ddddbe69c45d3c5c833d1"


def test_zenodo_doi_in_all_citation_surfaces():
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    workspace = (ROOT / "web" / "workspace.html").read_text(encoding="utf-8")
    assert f'doi: "{DOI}"' in cff
    assert 'version: "2.0.0"' in cff
    assert 'family-names: "Asamoah"' in cff
    assert 'given-names: "Benjamin O."' in cff
    assert f"https://doi.org/{DOI}" in readme
    assert f"https://doi.org/{DOI}" in workspace
    assert RELEASE_SHA in readme and RELEASE_SHA in workspace
    assert "id=\"citationText\"" in workspace
    citation_region = workspace.split('id="citationDialog"', 1)[1].split("</dialog>", 1)[0]
    assert f"https://doi.org/{DOI}" in citation_region
    assert "no DOI or peer-reviewed publication claimed" not in citation_region
