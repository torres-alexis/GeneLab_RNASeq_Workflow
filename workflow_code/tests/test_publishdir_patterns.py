"""Modules publish via entry-workflow output {}, not process publishDir."""

import re
from pathlib import Path

MODULES = Path(__file__).resolve().parents[1] / "modules"


def _hits(rx):
    out = []
    for path in sorted(MODULES.glob("*.nf")):
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if rx.search(line):
                out.append(f"{path.name}:{i}: {line.strip()}")
    return out


def test_no_publishdir():
    assert not _hits(re.compile(r"publishDir"))


def test_no_glob_string_pattern_closures():
    assert not _hits(re.compile(r"pattern:\s*\{\s*[\"']"))
