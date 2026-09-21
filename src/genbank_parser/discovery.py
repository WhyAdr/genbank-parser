"""Deterministic, symlink-safe discovery and raw-file fingerprinting."""

from __future__ import annotations

import hashlib
import os
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

GENBANK_SUFFIXES = (".gb", ".gbk", ".gbff")
GENBANK_GZIP_SUFFIXES = tuple(f"{suffix}.gz" for suffix in GENBANK_SUFFIXES)
DISCOVERABLE_SUFFIXES = frozenset((*GENBANK_SUFFIXES, *GENBANK_GZIP_SUFFIXES))


@dataclass(frozen=True)
class DiscoveredInput:
    """One source with both user-facing and canonical filesystem identity."""

    requested_path: str
    source_identity: str
    resolved_path: Path
    display_path: str
    discovery_root: Path | None


@dataclass(frozen=True)
class FileFingerprint:
    sha256: str
    size_bytes: int
    mtime_ns: int


def _canonical(path: Path) -> Path:
    return path.resolve(strict=True)


def _identity(path: Path) -> str:
    """Return the stable local identity used for persistent matching."""

    return os.path.normcase(os.fspath(path))


def _display(path: Path) -> str:
    return os.path.normpath(os.fspath(path))


def _is_excluded(path: Path, excluded: tuple[Path, ...]) -> bool:
    candidate = path.resolve(strict=False)
    for item in excluded:
        target = item.resolve(strict=False)
        if candidate == target:
            return True
        try:
            candidate.relative_to(target)
        except ValueError:
            continue
        return True
    return False


def _recognized(path: Path) -> bool:
    return path.name.casefold().endswith(tuple(DISCOVERABLE_SUFFIXES))


def discover_inputs(
    inputs: Iterable[str | Path],
    *,
    recursive: bool = True,
    exclude: Iterable[str | Path] = (),
) -> list[DiscoveredInput]:
    """Discover GenBank files in stable order.

    Explicit files are accepted even when their suffix is unusual because the
    canonical parser can identify gzip by magic bytes. Directory traversal is
    intentionally restricted to recognized GenBank suffixes.
    """

    excluded = tuple(Path(item) for item in exclude if item is not None)
    observations: dict[str, list[tuple[str, str, Path | None]]] = {}

    def observe(requested: str, path: Path, discovery_root: Path | None) -> None:
        resolved = _canonical(path)
        key = os.path.normcase(os.fspath(resolved))
        observations.setdefault(key, []).append(
            (requested, _display(path), discovery_root)
        )

    for raw_item in inputs:
        requested = os.fspath(raw_item)
        path = Path(requested)
        if path.is_file():
            if _is_excluded(path, excluded):
                continue
            observe(requested, path, None)
            continue
        if not path.is_dir():
            continue
        root = path.resolve(strict=True)
        candidates = root.rglob("*") if recursive else root.glob("*")
        for candidate in candidates:
            if not candidate.is_file() or not _recognized(candidate):
                continue
            if _is_excluded(candidate, excluded):
                continue
            observe(requested, candidate, root)
    seen: list[DiscoveredInput] = []
    for observations_for_source in observations.values():
        # Overlapping roots and explicit-file arguments can observe the same
        # source through different spellings.  Select metadata by a stable
        # total order rather than first-seen order so all argument
        # permutations produce the same logical source description.
        requested, display, discovery_root = min(
            observations_for_source,
            key=lambda item: (
                os.path.normcase(item[1]),
                item[1],
                os.path.normcase(item[0]),
                item[0],
                "" if item[2] is None else os.path.normcase(os.fspath(item[2])),
                "" if item[2] is None else os.fspath(item[2]),
            ),
        )
        resolved = _canonical(Path(display))
        seen.append(
            DiscoveredInput(
                requested_path=requested,
                source_identity=_identity(resolved),
                resolved_path=resolved,
                display_path=display,
                discovery_root=discovery_root,
            )
        )
    return sorted(
        seen,
        key=lambda item: (os.path.normcase(item.display_path), item.display_path),
    )


def fingerprint_file(path: str | Path, *, chunk_size: int = 1 << 20) -> FileFingerprint:
    """Return raw-byte SHA-256, size, and nanosecond mtime for a file."""

    target = Path(path)
    digest = hashlib.sha256()
    size = 0
    with target.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
            size += len(chunk)
    stat = target.stat()
    return FileFingerprint(digest.hexdigest(), size, stat.st_mtime_ns)


def snapshot_file(path: str | Path) -> tuple[bytes, FileFingerprint]:
    """Read a source once and return the exact bytes plus their fingerprint.

    Callers that parse the returned bytes can therefore bind provenance to the
    bytes that were actually analyzed instead of hashing a path and reopening
    it later.  The mtime is retained as contextual filesystem metadata; the
    digest and size describe the immutable snapshot.
    """

    target = Path(path)
    raw = target.read_bytes()
    stat = target.stat()
    return raw, FileFingerprint(
        sha256=hashlib.sha256(raw).hexdigest(),
        size_bytes=len(raw),
        mtime_ns=stat.st_mtime_ns,
    )


def _stem(path: Path) -> str:
    name = path.name
    folded = name.casefold()
    for suffix in sorted(DISCOVERABLE_SUFFIXES, key=len, reverse=True):
        if folded.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _safe_key(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    return cleaned or "sample"


def assign_sample_keys(
    discovered: Iterable[DiscoveredInput],
    *,
    sha256_by_identity: Mapping[str, str] | None = None,
) -> dict[DiscoveredInput, str]:
    """Assign readable, deterministic, collision-safe cohort sample keys."""

    items = tuple(
        sorted(
            discovered,
            key=lambda item: (
                os.path.normcase(item.source_identity),
                item.source_identity,
                os.path.normcase(item.display_path),
                item.display_path,
            ),
        )
    )

    def digest_for(item: DiscoveredInput) -> str:
        if sha256_by_identity is not None:
            digest = sha256_by_identity.get(item.source_identity)
            if digest:
                return digest
        return fingerprint_file(item.resolved_path).sha256

    base_names = [_safe_key(_stem(item.resolved_path)) for item in items]
    counts: dict[str, int] = {}
    for base in base_names:
        counts[base.casefold()] = counts.get(base.casefold(), 0) + 1
    result: dict[DiscoveredInput, str] = {}
    used_keys: set[str] = set()
    for item, base in zip(items, base_names, strict=True):
        digest = ""
        if counts[base.casefold()] == 1:
            candidate = base
        else:
            digest = digest_for(item)[:8]
            candidate = f"{base}-{digest}"
        suffix = 2
        while candidate.casefold() in used_keys:
            if not digest:
                digest = digest_for(item)[:8]
            candidate = f"{base}-{digest}-{suffix}"
            suffix += 1
        used_keys.add(candidate.casefold())
        result[item] = candidate
    return result


__all__ = [
    "DISCOVERABLE_SUFFIXES",
    "GENBANK_GZIP_SUFFIXES",
    "GENBANK_SUFFIXES",
    "DiscoveredInput",
    "FileFingerprint",
    "assign_sample_keys",
    "discover_inputs",
    "fingerprint_file",
    "snapshot_file",
]
