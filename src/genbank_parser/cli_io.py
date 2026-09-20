"""Small, shared CLI I/O and publication helpers.

The command modules deliberately keep their analysis code independent from
these helpers.  The unified CLI uses this module to make output publication
predictable: data are written atomically, diagnostics stay on stderr, and an
input file can never be selected as its own output.
"""

from __future__ import annotations

import contextlib
import gzip
import os
import shutil
import sys
import tempfile
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


class GbparseError(Exception):
    """Base class for expected command failures."""


class InputError(GbparseError):
    """An input could not be read or parsed."""


class OutputError(GbparseError):
    """An output could not be safely published."""


class SerializationError(OutputError):
    """A renderer produced an invalid or unpublishable result."""


@dataclass(frozen=True)
class InputSource:
    """An opened source with a stable human-readable label."""

    label: str
    path: Path | None
    stream: TextIO
    compressed: bool = False

    @property
    def source_kind(self) -> str:
        if self.path is not None:
            return "path"
        return "stdin" if self.label == "-" else "stream"


def eprint(message: object) -> None:
    """Write a diagnostic to stderr."""

    print(message, file=sys.stderr)


def _normal_path(path: str | Path) -> str:
    # normcase matters on Windows, where two differently cased spellings can
    # still designate the same input/output file.
    # Resolve existing parent links too, so an output through a symlinked
    # directory cannot alias a path-backed input that has not been created yet.
    try:
        resolved = Path(path).resolve(strict=False)
    except OSError:
        resolved = Path(os.path.abspath(os.fspath(path)))
    return os.path.normcase(os.fspath(resolved))


def paths_same(left: str | Path, right: str | Path) -> bool:
    """Return whether two paths resolve to the same filesystem object."""

    left_path = Path(left)
    right_path = Path(right)
    try:
        if left_path.exists() and right_path.exists():
            return os.path.samefile(left_path, right_path)
    except OSError:
        # The normalized absolute comparison below remains useful for a
        # not-yet-created output or a filesystem without samefile support.
        pass
    return _normal_path(left_path) == _normal_path(right_path)


def reject_input_output_collision(
    inputs: Iterable[str | Path | None], output: str | Path | None
) -> None:
    """Reject an output path that aliases any input path.

    ``resolve(strict=False)`` and ``samefile`` together cover ordinary path
    aliases and an existing symlink to an input without requiring the output
    to exist already.
    """

    if output is None or os.fspath(output) == "-":
        return
    output_path = Path(output)
    for source in inputs:
        if source is None or os.fspath(source) == "-":
            continue
        if paths_same(source, output_path):
            raise OutputError(
                f"output path must not overwrite input: {output_path}"
            )


@contextlib.contextmanager
def open_input(
    source: str | Path | TextIO,
    *,
    stdin: TextIO | None = None,
) -> Iterator[InputSource]:
    """Open a text source for callers that need an explicit source object.

    Parsing itself lives in :func:`genbank_parser.io.read_genbank`; this
    helper is intentionally small and is mainly useful to exporters or CLI
    code that needs a source label.  File-like objects are not closed by this
    context manager.
    """

    if hasattr(source, "read"):
        stream = source  # type: ignore[assignment]
        label = str(getattr(stream, "name", "<stream>"))
        yield InputSource(label=label, path=None, stream=stream)
        return

    source_text = os.fspath(source)
    if source_text == "-":
        stream = stdin or sys.stdin
        yield InputSource(label="-", path=None, stream=stream)
        return

    path = Path(source)
    try:
        with path.open("rb") as probe:
            compressed = probe.read(2) == b"\x1f\x8b"
        handle = (
            gzip.open(path, "rt", encoding="utf-8", errors="replace")  # noqa: SIM115
            if compressed
            else path.open("r", encoding="utf-8", errors="replace")
        )
    except OSError as exc:
        raise InputError(f"could not open input {path}: {exc}") from exc
    try:
        yield InputSource(
            label=str(path),
            path=path,
            stream=handle,
            compressed=compressed,
        )
    finally:
        handle.close()


@contextlib.contextmanager
def atomic_text_writer(
    path: str | Path,
    *,
    force: bool = False,
    inputs: Iterable[str | Path | None] = (),
) -> Iterator[TextIO]:
    """Yield a UTF-8 text handle and atomically replace ``path`` on success."""

    output = Path(path)
    reject_input_output_collision(inputs, output)
    if output.exists() and not force:
        raise OutputError(f"output already exists; use --force: {output}")
    if output.exists() and output.is_dir():
        raise OutputError(f"output path is a directory: {output}")
    temporary_name: str | None = None
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            delete=False,
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
        ) as handle:
            temporary_name = handle.name
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, output)
        temporary_name = None
    except OSError as exc:
        raise OutputError(f"could not publish output {output}: {exc}") from exc
    finally:
        if temporary_name:
            try:
                Path(temporary_name).unlink(missing_ok=True)
            except OSError:
                pass


def write_text(
    rendered: str,
    output: str | Path | None = None,
    *,
    force: bool = False,
    inputs: Iterable[str | Path | None] = (),
) -> None:
    """Write rendered data to stdout or atomically to a file."""

    if output is None or os.fspath(output) == "-":
        sys.stdout.write(rendered)
        return
    with atomic_text_writer(output, force=force, inputs=inputs) as handle:
        handle.write(rendered)


@contextlib.contextmanager
def atomic_bytes_writer(
    path: str | Path,
    *,
    force: bool = False,
    inputs: Iterable[str | Path | None] = (),
) -> Iterator[object]:
    """Yield a binary handle and atomically replace ``path`` on success."""

    output = Path(path)
    reject_input_output_collision(inputs, output)
    if output.exists() and not force:
        raise OutputError(f"output already exists; use --force: {output}")
    if output.exists() and output.is_dir():
        raise OutputError(f"output path is a directory: {output}")
    temporary_name: str | None = None
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="wb",
            delete=False,
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
        ) as handle:
            temporary_name = handle.name
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, output)
        temporary_name = None
    except OSError as exc:
        raise OutputError(f"could not publish output {output}: {exc}") from exc
    finally:
        if temporary_name:
            try:
                Path(temporary_name).unlink(missing_ok=True)
            except OSError:
                pass


def write_bytes(
    rendered: bytes,
    output: str | Path,
    *,
    force: bool = False,
    inputs: Iterable[str | Path | None] = (),
) -> None:
    """Write bytes atomically to a file."""

    with atomic_bytes_writer(output, force=force, inputs=inputs) as handle:
        handle.write(rendered)  # type: ignore[union-attr]


def _write_staged_file(path: Path, payload: str | bytes) -> None:
    mode = "wb" if isinstance(payload, bytes) else "w"
    kwargs = {"encoding": "utf-8", "newline": ""} if mode == "w" else {}
    with path.open(mode, **kwargs) as handle:  # type: ignore[arg-type]
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def publish_staged_files(
    replacements: Sequence[tuple[str | Path, str | Path]],
    *,
    force: bool = False,
    inputs: Iterable[str | Path | None] = (),
) -> None:
    """Install several already-staged files as one rollback-protected set.

    All destinations are preflighted before the first replacement. Existing
    destinations are moved to sibling backups and restored if any later
    replacement fails, so callers never observe only a prefix of a logical
    multi-file publication.
    """

    items = tuple((Path(staged), Path(destination)) for staged, destination in replacements)
    if not items:
        return
    destinations = tuple(destination for _staged, destination in items)
    for index, destination in enumerate(destinations):
        if any(paths_same(destination, other) for other in destinations[index + 1 :]):
            raise OutputError(f"publication destinations must be distinct: {destination}")
        reject_input_output_collision(inputs, destination)
        if destination.exists() and destination.is_dir():
            raise OutputError(f"output path is a directory: {destination}")
    for staged, _destination in items:
        if not staged.is_file():
            raise OutputError(f"staged output is missing: {staged}")
    if not force:
        existing = next((destination for destination in destinations if destination.exists()), None)
        if existing is not None:
            raise OutputError(f"output already exists; use --force: {existing}")

    backups: dict[Path, Path] = {}
    installed: list[Path] = []
    try:
        for _staged, destination in items:
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists():
                handle, backup_name = tempfile.mkstemp(
                    prefix=f".{destination.name}.",
                    suffix=".backup",
                    dir=destination.parent,
                )
                os.close(handle)
                backup = Path(backup_name)
                backup.unlink(missing_ok=True)
                os.replace(destination, backup)
                backups[destination] = backup
        for staged, destination in items:
            os.replace(staged, destination)
            installed.append(destination)
    except OSError as exc:
        for destination in installed:
            try:
                destination.unlink(missing_ok=True)
            except OSError:
                pass
        for destination, backup in backups.items():
            if backup.exists() and not destination.exists():
                try:
                    os.replace(backup, destination)
                except OSError:
                    pass
        raise OutputError(f"could not publish output set: {exc}") from exc
    finally:
        for staged, _destination in items:
            staged.unlink(missing_ok=True)
        for backup in backups.values():
            backup.unlink(missing_ok=True)


def publish_file_set(
    files: Mapping[str | Path, str | bytes],
    *,
    force: bool = False,
    inputs: Iterable[str | Path | None] = (),
) -> None:
    """Render and publish a set of files through :func:`publish_staged_files`."""

    if not files:
        return
    staged: list[tuple[Path, Path]] = []
    try:
        for raw_destination, payload in files.items():
            destination = Path(raw_destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            handle, staged_name = tempfile.mkstemp(
                prefix=f".{destination.name}.",
                suffix=".tmp",
                dir=destination.parent,
            )
            os.close(handle)
            staged_path = Path(staged_name)
            _write_staged_file(staged_path, payload)
            staged.append((staged_path, destination))
        publish_staged_files(staged, force=force, inputs=inputs)
        staged = []
    finally:
        for staged_path, _destination in staged:
            staged_path.unlink(missing_ok=True)


def publish_directory(
    files: dict[str, str | bytes],
    output_dir: str | Path,
    *,
    force: bool = False,
    inputs: Iterable[str | Path | None] = (),
) -> Path:
    """Atomically publish a small multi-file export directory.

    Files are assembled in a sibling staging directory.  When ``force`` is
    used for an existing destination, the old directory is moved aside until
    the new directory has been installed, allowing a failed publication to
    restore it.
    """

    destination = Path(output_dir)
    reject_input_output_collision(inputs, destination)
    if destination.exists() and destination.is_dir():
        destination_resolved = destination.resolve()
        for source in inputs:
            if source is None or os.fspath(source) == "-":
                continue
            source_path = Path(source)
            try:
                source_path.resolve().relative_to(destination_resolved)
            except ValueError:
                continue
            raise OutputError(
                f"output directory would replace input: {destination}"
            )
    if destination.exists() and not destination.is_dir():
        raise OutputError(f"output directory is not a directory: {destination}")
    if destination.exists() and not force:
        raise OutputError(f"output directory already exists; use --force: {destination}")
    staging: Path | None = None
    backup: Path | None = None
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(
            tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent)
        )
        for relative_name, payload in sorted(files.items()):
            relative = Path(relative_name)
            if relative.is_absolute() or ".." in relative.parts:
                raise OutputError(f"invalid export filename: {relative_name}")
            target = staging / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            _write_staged_file(target, payload)

        if destination.exists():
            backup = Path(
                tempfile.mkdtemp(prefix=f".{destination.name}.backup.", dir=destination.parent)
            )
            backup.rmdir()
            os.replace(destination, backup)
        try:
            assert staging is not None
            os.replace(staging, destination)
        except OSError:
            if backup is not None and not destination.exists():
                os.replace(backup, destination)
                backup = None
            raise
        staging = None
        if backup is not None:
            shutil.rmtree(backup)
            backup = None
        return destination
    except OSError as exc:
        raise OutputError(f"could not publish output directory {destination}: {exc}") from exc
    finally:
        if staging is not None and staging.exists():
            shutil.rmtree(staging, ignore_errors=True)
        if backup is not None and backup.exists() and not destination.exists():
            try:
                os.replace(backup, destination)
            except OSError:
                pass


def publish_directory_tree(
    staging_dir: str | Path,
    output_dir: str | Path,
    *,
    force: bool = False,
    inputs: Iterable[str | Path | None] = (),
) -> Path:
    """Atomically install an already assembled directory tree."""

    staging = Path(staging_dir)
    destination = Path(output_dir)
    if not staging.is_dir():
        raise OutputError(f"staging directory does not exist: {staging}")
    reject_input_output_collision(inputs, destination)
    if destination.exists() and not destination.is_dir():
        raise OutputError(f"output directory is not a directory: {destination}")
    if destination.exists() and not force:
        raise OutputError(f"output directory already exists; use --force: {destination}")
    backup: Path | None = None
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            backup = Path(
                tempfile.mkdtemp(prefix=f".{destination.name}.backup.", dir=destination.parent)
            )
            backup.rmdir()
            os.replace(destination, backup)
        os.replace(staging, destination)
        if backup is not None:
            shutil.rmtree(backup)
        return destination
    except OSError as exc:
        if backup is not None and backup.exists() and not destination.exists():
            try:
                os.replace(backup, destination)
                backup = None
            except OSError:
                pass
        raise OutputError(f"could not publish output directory {destination}: {exc}") from exc
    finally:
        if backup is not None and backup.exists() and not destination.exists():
            shutil.rmtree(backup, ignore_errors=True)


def json_text(value: object) -> str:
    """Serialize JSON with stable ordering and a final newline."""

    try:
        return (
            __import__("json").dumps(
                value,
                ensure_ascii=False,
                allow_nan=False,
                sort_keys=True,
                indent=2,
            )
            + "\n"
        )
    except (TypeError, ValueError) as exc:
        raise SerializationError(f"could not serialize JSON: {exc}") from exc


__all__ = [
    "GbparseError",
    "InputError",
    "InputSource",
    "OutputError",
    "SerializationError",
    "atomic_bytes_writer",
    "atomic_text_writer",
    "eprint",
    "json_text",
    "open_input",
    "paths_same",
    "publish_directory",
    "publish_directory_tree",
    "publish_file_set",
    "publish_staged_files",
    "reject_input_output_collision",
    "write_bytes",
    "write_text",
]
