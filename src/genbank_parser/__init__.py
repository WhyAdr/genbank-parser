"""GenBank Feature Parser & Genome Annotation Query Engine."""

from __future__ import annotations

from .io import (
    extract_xref_sources,
    extract_xrefs,
    get_notes,
    get_qual,
    iter_genbank,
    parse_features,
    read_genbank,
)
from .model import GenBankDocument, GenBankFeature, GenBankRecord
from .neighborhood import NeighborhoodResult, build_neighborhood
from .operons import OperonCluster, OperonPair, OperonResult, build_operon_result

# Keep the source-checkout CLI and an installed wheel on the same release
# value.  The project metadata in pyproject.toml is intentionally static too.
__version__ = "0.9.2"

__all__ = [
    "GenBankDocument",
    "GenBankFeature",
    "GenBankRecord",
    "NeighborhoodResult",
    "OperonCluster",
    "OperonPair",
    "OperonResult",
    "__version__",
    "build_neighborhood",
    "build_operon_result",
    "extract_xref_sources",
    "extract_xrefs",
    "get_notes",
    "get_qual",
    "iter_genbank",
    "parse_features",
    "read_genbank",
]
