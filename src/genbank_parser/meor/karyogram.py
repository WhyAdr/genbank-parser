"""Record-length-aware ASCII density karyograms."""
from __future__ import annotations

from typing import Any

from ..model import GenBankRecord
from .models import MeorHit


def generate_meor_karyograms(
    records: list[GenBankRecord] | tuple[GenBankRecord, ...],
    hits: list[MeorHit] | tuple[MeorHit, ...],
    *,
    window_size: int = 50_000,
) -> dict[str, dict[str, Any]]:
    if window_size <= 0:
        raise ValueError("window_size must be positive")
    karyograms: dict[str, dict[str, Any]] = {}
    for record in records:
        num_windows = max(1, (record.length + window_size - 1) // window_size)
        counts = [0] * num_windows
        record_hits = [hit for hit in hits if hit.contig == record.id]
        for hit in record_hits:
            index = min(num_windows - 1, (hit.start - 1) // window_size)
            counts[index] += 1
        karyograms[record.id] = {
            "span_bp": record.length,
            "window_size": window_size,
            "total_hits": len(record_hits),
            "ascii_map": "".join("#" if count else "-" for count in counts),
        }
    return karyograms
