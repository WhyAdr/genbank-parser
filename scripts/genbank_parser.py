#!/usr/bin/env python3
import sys
from pathlib import Path

# Prevent scripts/ directory from shadowing the genbank_parser package
_scripts_dir = str(Path(__file__).resolve().parent)
while _scripts_dir in sys.path:
    sys.path.remove(_scripts_dir)

_src = str(Path(__file__).resolve().parent.parent / "src")
if _src not in sys.path:
    sys.path.insert(0, _src)

from genbank_parser.io import extract_xrefs, get_notes, get_qual, parse_features, read_genbank
from genbank_parser.model import GenBankDocument, GenBankFeature, GenBankRecord

__all__ = [
    "read_genbank",
    "parse_features",
    "get_qual",
    "get_notes",
    "extract_xrefs",
    "GenBankDocument",
    "GenBankRecord",
    "GenBankFeature",
]
