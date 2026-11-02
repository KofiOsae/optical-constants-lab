"""Static regression check for filename-derived thickness in both interfaces."""
from pathlib import Path
import re
import pytest

@pytest.mark.parametrize("page",["workspace.html","multilayer.html"])
def test_filename_regex_is_not_overescaped(page):
    source=(Path(__file__).resolve().parents[1]/"web"/page).read_text()
    match=re.search(r"function oclThicknessFromName\(name\)\{([^\n]*)",source)
    assert match, "Thickness function missing"
    body=match.group(1)
    assert r"\s*" in body and r"\." in body
    assert r"\\s" not in body and r"\\." not in body
