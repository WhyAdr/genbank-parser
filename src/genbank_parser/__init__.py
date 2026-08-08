"""GenBank Feature Parser & Genome Annotation Query Engine."""
from __future__ import annotations

from .io import extract_xrefs, get_notes, get_qual, parse_features, read_genbank
from .model import GenBankDocument, GenBankFeature, GenBankRecord

__version__ = "0.2.0"

__all__ = [
    "read_genbank",
    "parse_features",
    "get_qual",
    "get_notes",
    "extract_xrefs",
    "GenBankDocument",
    "GenBankRecord",
    "GenBankFeature",
    "__version__",
]
