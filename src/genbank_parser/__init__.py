"""GenBank Feature Parser & Genome Annotation Query Engine."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from .io import extract_xrefs, get_notes, get_qual, parse_features, read_genbank
from .model import GenBankDocument, GenBankFeature, GenBankRecord

try:
    __version__ = version("genbank-parser")
except PackageNotFoundError:
    # Source checkouts can import the package before installation; the version
    # remains defined solely by pyproject metadata in installed environments.
    __version__ = "unknown"

__all__ = [
    "GenBankDocument",
    "GenBankFeature",
    "GenBankRecord",
    "__version__",
    "extract_xrefs",
    "get_notes",
    "get_qual",
    "parse_features",
    "read_genbank",
]
