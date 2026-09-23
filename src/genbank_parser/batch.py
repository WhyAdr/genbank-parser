"""Manifest-driven execution of registered single-input gbparse commands."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import secrets
import shutil
import subprocess
import sys
import tempfile
import time
from argparse import Namespace
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from importlib import resources
from pathlib import Path
from typing import Literal

from . import __version__
from .cli_io import (
    InputError,
    OutputError,
    OutputRecoveryError,
    publish_directory_tree,
    reject_input_output_collision,
    write_text,
)
from .discovery import (
    DiscoveredInput,
    assign_sample_keys,
    discover_inputs,
    snapshot_file,
)

OutputMode = Literal["file", "directory", "multi_file"]


class BatchUsageError(ValueError):
    """The requested batch command/adapter combination is invalid."""


@dataclass(frozen=True)
class ExpectedOutput:
    relative_path: str


@dataclass(frozen=True)
class BatchCommandSpec:
    name: str
    input_arity: Literal[1] = 1
    output_mode: OutputMode = "file"
    default_suffix: str | None = ".out"
    batchable: bool = True
    forbidden_options: tuple[str, ...] = ("--output", "--output-dir", "--force")
    build_job_argv: Callable[[str, str, Sequence[str], Namespace], list[str]] = field(
        default=lambda command, input_path, tail, _options: [command, input_path, *tail]
    )
    collect_outputs: Callable[[Path, str, Sequence[str], Namespace], tuple[ExpectedOutput, ...]] = field(
        default=lambda sample_dir, input_path, tail, _options: (ExpectedOutput("result.out"),)
    )
    capture_extra_outputs: bool = False


@dataclass(frozen=True)
class BatchExecutionResult:
    manifest: dict[str, object]
    failed_count: int
    exit_code: int


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _sha256(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def _format_suffix(tail: Sequence[str], fallback: str, options: Namespace | None = None) -> str:
    emit = getattr(options, "emit", None) if options is not None else None
    if emit in {"faa", "ffn"}:
        return "." + emit
    value = getattr(options, "format", None) if options is not None else None
    if options is not None and getattr(options, "json", False):
        value = "json"
    if value is None:
        for index, token in enumerate(tail):
            if token == "--format" and index + 1 < len(tail):
                value = tail[index + 1]
                break
            if token.startswith("--format="):
                value = token.split("=", 1)[1]
                break
        if value is None and "--json" in tail:
            value = "json"
    if value is None:
        return fallback
    return {
        "json": ".json",
        "jsonl": ".jsonl",
        "tsv": ".tsv",
        "csv": ".csv",
        "text": ".txt",
        "annotations-tsv": ".tsv",
        "faa": ".faa",
        "ffn": ".ffn",
        "fna": ".fna",
        "gff3": ".gff3",
        "bed12": ".bed",
    }.get(value, fallback)


def _option_value(tail: Sequence[str], name: str) -> str | None:
    for index, token in enumerate(tail):
        if token == name and index + 1 < len(tail):
            return tail[index + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def _source_stem(input_path: str) -> str:
    name = Path(input_path).name
    if name.casefold().endswith(".gz"):
        name = Path(name).stem
    return Path(name).stem


def _file_spec(name: str, suffix: str) -> BatchCommandSpec:
    def build(command: str, input_path: str, tail: Sequence[str], options: Namespace) -> list[str]:
        output_name = f"result{_format_suffix(tail, suffix, options)}"
        return [command, input_path, *tail, "--output", "{OUTPUT}/" + output_name]

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str], options: Namespace) -> tuple[ExpectedOutput, ...]:
        return (ExpectedOutput(f"result{_format_suffix(tail, suffix, options)}"),)

    return BatchCommandSpec(name=name, default_suffix=suffix, build_job_argv=build, collect_outputs=collect)


def _positional_spec(name: str, suffix: str) -> BatchCommandSpec:
    def build(command: str, input_path: str, tail: Sequence[str], options: Namespace) -> list[str]:
        return [command, input_path, *tail, "{OUTPUT}/result" + _format_suffix(tail, suffix, options)]

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str], options: Namespace) -> tuple[ExpectedOutput, ...]:
        return (ExpectedOutput("result" + _format_suffix(tail, suffix, options)),)

    return BatchCommandSpec(
        name=name,
        default_suffix=suffix,
        build_job_argv=build,
        collect_outputs=collect,
        forbidden_options=("--output", "--output-dir", "--force"),
    )


def _sequence_spec() -> BatchCommandSpec:
    def build(command: str, input_path: str, tail: Sequence[str], options: Namespace) -> list[str]:
        return [command, input_path, *tail, "--fna", "{OUTPUT}/genome.fna", "--ffn", "{OUTPUT}/cds.ffn"]

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str], options: Namespace) -> tuple[ExpectedOutput, ...]:
        return (ExpectedOutput("genome.fna"), ExpectedOutput("cds.ffn"))

    return BatchCommandSpec(
        name="sequence",
        output_mode="multi_file",
        default_suffix=None,
        forbidden_options=("--output", "--output-dir", "--force", "--fna", "--ffn"),
        build_job_argv=build,
        collect_outputs=collect,
    )


def _neighborhood_spec() -> BatchCommandSpec:
    def build(command: str, input_path: str, tail: Sequence[str], options: Namespace) -> list[str]:
        argv = [command, input_path, *tail, "--output", "{OUTPUT}/result" + _format_suffix(tail, ".txt", options)]
        if "--visualize" in tail or "--viz-output" in tail:
            argv.extend(("--visualize", "--viz-output", "{OUTPUT}/neighborhood.svg"))
        return argv

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str], options: Namespace) -> tuple[ExpectedOutput, ...]:
        result = [ExpectedOutput("result" + _format_suffix(tail, ".txt", options))]
        if "--visualize" in tail or "--viz-output" in tail:
            result.append(ExpectedOutput("neighborhood.svg"))
        return tuple(result)

    return BatchCommandSpec(
        name="neighborhood",
        output_mode="multi_file",
        default_suffix=None,
        forbidden_options=("--output", "--output-dir", "--force", "--viz-output"),
        build_job_argv=build,
        collect_outputs=collect,
    )


def _phylo_spec() -> BatchCommandSpec:
    def build(command: str, input_path: str, tail: Sequence[str], options: Namespace) -> list[str]:
        return [command, input_path, *tail, "--output", "{OUTPUT}/report" + _format_suffix(tail, ".txt", options), "--output-dir", "{OUTPUT}/markers"]

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str], options: Namespace) -> tuple[ExpectedOutput, ...]:
        return (ExpectedOutput("report" + _format_suffix(tail, ".txt", options)),)

    return BatchCommandSpec(
        name="phylo",
        output_mode="multi_file",
        default_suffix=None,
        capture_extra_outputs=True,
        forbidden_options=("--output", "--output-dir", "--force"),
        build_job_argv=build,
        collect_outputs=collect,
    )


def _export_spec() -> BatchCommandSpec:
    def build(command: str, input_path: str, tail: Sequence[str], options: Namespace) -> list[str]:
        if getattr(options, "format", None) == "ncbi-table" or _option_value(tail, "--format") == "ncbi-table":
            return [command, input_path, *tail, "--output-dir", "{OUTPUT}"]
        suffix = {"annotations-tsv": ".tsv", "jsonl": ".jsonl", "faa": ".faa", "ffn": ".ffn", "fna": ".fna", "gff3": ".gff3", "bed12": ".bed"}.get(getattr(options, "format", None), ".out")
        return [command, input_path, *tail, "--output", "{OUTPUT}/result" + suffix]

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str], options: Namespace) -> tuple[ExpectedOutput, ...]:
        if getattr(options, "format", None) == "ncbi-table" or _option_value(tail, "--format") == "ncbi-table":
            prefix = _source_stem(input_path)
            return tuple(
                ExpectedOutput(name)
                for name in (
                    f"{prefix}.fsa",
                    f"{prefix}.tbl",
                    f"{prefix}.export.json",
                )
            )
        suffix = {"annotations-tsv": ".tsv", "jsonl": ".jsonl", "faa": ".faa", "ffn": ".ffn", "fna": ".fna", "gff3": ".gff3", "bed12": ".bed"}.get(getattr(options, "format", None), ".out")
        return (ExpectedOutput("result" + suffix),)

    return BatchCommandSpec(
        name="export",
        output_mode="file",
        default_suffix=None,
        forbidden_options=("--output", "--output-dir", "--force"),
        build_job_argv=build,
        collect_outputs=collect,
    )


_REGISTRY: dict[str, BatchCommandSpec] = {
    "validate": _file_spec("validate", ".txt"),
    "summary": _file_spec("summary", ".txt"),
    "metadata": _file_spec("metadata", ".txt"),
    "extract": _positional_spec("extract", ".tsv"),
    "search": _file_spec("search", ".txt"),
    "locus": _file_spec("locus", ".txt"),
    "neighborhood": _neighborhood_spec(),
    "region": _file_spec("region", ".gbff"),
    "fasta": _positional_spec("fasta", ".faa"),
    "sequence": _sequence_spec(),
    "codon": _file_spec("codon", ".tsv"),
    "functional": _file_spec("functional", ".tsv"),
    "discover": _file_spec("discover", ".txt"),
    "phylo": _phylo_spec(),
    "crispr": _file_spec("crispr", ".txt"),
    "gff": _positional_spec("gff", ".gff3"),
    "meor": _file_spec("meor", ".txt"),
    "mobilome": _file_spec("mobilome", ".txt"),
    "query": _file_spec("query", ".txt"),
    "operons": _file_spec("operons", ".txt"),
    "export": _export_spec(),
}

_UNSUPPORTED = {
    "compare": "compare takes multiple genomes; invoke it directly for a cohort matrix",
    "diff": "diff takes two annotation versions; invoke it directly",
    "batch-summary": "batch-summary is already a cohort aggregation; invoke it directly",
    "batch": "nested batch execution is not supported",
    "index": "index lifecycle commands are not single-input jobs",
    "records": "records split manages its own directory; use records directly",
}

_NON_BATCHABLE_SPECS = {
    name: BatchCommandSpec(
        name=name,
        batchable=False,
        output_mode="file",
        default_suffix=None,
    )
    for name in _UNSUPPORTED
}

# Custom auxiliary resources are intentionally not batchable until their raw
# bytes can be fingerprinted, snapshotted, published, and validated on resume.
# Keep this contract local to batch so direct commands retain their existing
# support for reviewed or user-supplied resources.
_EXTERNAL_RESOURCE_FLAGS: dict[str, tuple[str, ...]] = {
    "discover": ("--rules",),
    "meor": ("--markers", "--pathways"),
    "mobilome": ("--database-dir",),
}


def registered_commands() -> dict[str, BatchCommandSpec]:
    return {**_REGISTRY, **_NON_BATCHABLE_SPECS}


def get_command_spec(command: str, tail: Sequence[str] = ()) -> BatchCommandSpec:
    if command in _UNSUPPORTED:
        raise BatchUsageError(f"unsupported batch command {command!r}: {_UNSUPPORTED[command]}")
    try:
        spec = registered_commands()[command]
    except KeyError as exc:
        raise BatchUsageError(f"unsupported batch command {command!r}; use a registered single-input command") from exc
    if not spec.batchable:
        raise BatchUsageError(f"unsupported batch command {command!r}: {_UNSUPPORTED[command]}")
    if command == "export" and _option_value(tail, "--format") == "ncbi-table":
        return replace(spec, output_mode="directory")
    return spec


def _validate_tail(spec: BatchCommandSpec, tail: Sequence[str]) -> None:
    for token in tail:
        for forbidden in spec.forbidden_options:
            if token == forbidden or token.startswith(forbidden + "="):
                raise BatchUsageError(
                    f"batch owns {forbidden}; remove it from the wrapped command tail"
                )


def _option_present(tail: Sequence[str], flag: str) -> bool:
    return any(token == flag or token.startswith(flag + "=") for token in tail)


def _validate_batch_dependencies(command: str, tail: Sequence[str]) -> None:
    """Reject auxiliary resources whose provenance is not batch-bound yet."""

    for flag in _EXTERNAL_RESOURCE_FLAGS.get(command, ()):
        if _option_present(tail, flag):
            raise BatchUsageError(
                f"batch command {command!r} does not accept {flag}; "
                "run the command directly until auxiliary-resource snapshot "
                "provenance is implemented"
            )


def _validate_command_argv(command: str, argv: Sequence[str]) -> Namespace:
    from .cli import create_command_parser

    parser = create_command_parser(command)
    try:
        return parser.parse_args(list(argv))
    except SystemExit as exc:
        raise BatchUsageError(
            f"batch command {command!r} arguments failed pre-validation"
        ) from exc


def _build_argv(
    spec: BatchCommandSpec,
    command: str,
    input_path: Path,
    sample_dir: Path,
    tail: Sequence[str],
    options: Namespace,
) -> list[str]:
    raw = spec.build_job_argv(command, str(input_path), tail, options)
    return [
        sys.executable,
        "-m",
        "genbank_parser.cli",
        *[str(sample_dir) if value == "{OUTPUT}" else str(sample_dir / value[len("{OUTPUT}/") :]) if value.startswith("{OUTPUT}/") else value for value in raw],
    ]


def _logical_argv(
    spec: BatchCommandSpec,
    command: str,
    source: DiscoveredInput,
    sample_key: str,
    tail: Sequence[str],
    options: Namespace,
    final_output_dir: Path,
) -> list[str]:
    """Build replayable argv without physical snapshot/work-tree paths."""

    argv = _build_argv(
        spec,
        command,
        source.resolved_path.resolve(strict=False),
        final_output_dir.resolve(strict=False) / "outputs" / sample_key,
        tail,
        options,
    )
    # A published run already contains the recorded outputs.  Replay is a
    # regeneration command, so it must be safe to execute without first
    # deleting those outputs manually.
    return [*argv, "--force"]


def _scrub_diagnostics(
    value: object,
    *,
    stage_root: Path,
    snapshot_root: Path,
) -> str:
    """Remove disposable batch paths from persisted child diagnostics."""

    text = str(value)
    replacements = {
        str(stage_root): "<batch-output>",
        str(stage_root.resolve(strict=False)): "<batch-output>",
        str(snapshot_root): "<batch-snapshot>",
        str(snapshot_root.resolve(strict=False)): "<batch-snapshot>",
    }
    for raw, replacement in replacements.items():
        if not raw:
            continue
        for spelling in {raw, raw.replace("\\", "/"), raw.replace("/", "\\")}:
            text = text.replace(spelling, replacement)
    # Keep the persisted-log contract robust even when a platform or child
    # formats a path with a different absolute/relative prefix.
    return (
        text.replace(".gbparse-input.", "<batch-snapshot>.")
        .replace(".gbparse-inprogress", "<batch-output>")
        .replace(".work.", "<batch-work>.")
    )


def _manifest_path(run_dir: Path) -> Path:
    return run_dir / "batch-manifest.json"


def _write_manifest(path: Path, manifest: dict[str, object]) -> None:
    write_text(json.dumps(manifest, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2) + "\n", path, force=True)


def _job_id(source: DiscoveredInput) -> str:
    return "job-" + hashlib.sha256(source.source_identity.encode("utf-8")).hexdigest()[:16]


def _input_payload(
    source: DiscoveredInput,
    fingerprint=None,
    *,
    fallback: dict[str, object] | None = None,
) -> dict[str, object]:
    if fingerprint is None:
        try:
            fingerprint = snapshot_file(source.resolved_path)[1]
        except OSError:
            if fallback is None:
                raise
            payload = dict(fallback)
            payload.update(
                {
                    "requested_path": source.requested_path,
                    "source_identity": source.source_identity,
                    "display_path": source.display_path,
                    "resolved_path": str(source.resolved_path),
                }
            )
            return payload
    fp = fingerprint
    return {
        "requested_path": source.requested_path,
        "source_identity": source.source_identity,
        "display_path": source.display_path,
        "resolved_path": str(source.resolved_path),
        "sha256": fp.sha256,
        "size_bytes": fp.size_bytes,
        "mtime_ns": fp.mtime_ns,
    }


def _manifest_template(command: str, tail: Sequence[str], inputs: Sequence[DiscoveredInput], keys: dict[DiscoveredInput, str], jobs: int, on_error: str) -> dict[str, object]:
    return {
        "schema_version": "gbparse.batch.v1",
        "gbparse_version": __version__,
        "python_version": platform.python_version(),
        "biopython_version": __import__("Bio").__version__,
        "platform": platform.platform(),
        "created_at": _now(),
        "updated_at": _now(),
        "command": command,
        "command_args": list(tail),
        "runner": {
            "jobs": jobs,
            "on_error": on_error,
            "discovery_roots": sorted(
                {str(item.discovery_root) for item in inputs if item.discovery_root is not None}
            ),
            "working_directory": str(Path.cwd().resolve()),
        },
        "jobs": [
            {"job_id": _job_id(item), "sample_key": keys[item], "status": "pending", "input": _input_payload(item)}
            for item in inputs
        ],
    }


def _outputs_payload(sample_root: Path, expected: Sequence[ExpectedOutput]) -> list[dict[str, object]]:
    outputs: list[dict[str, object]] = []
    for item in expected:
        path = sample_root / item.relative_path
        if not path.is_file():
            raise OutputError(f"expected batch output is missing: {path}")
        digest, size = _sha256(path)
        relative_path = str(item.relative_path).replace("\\", "/")
        outputs.append({"path": str(path.relative_to(sample_root.parent.parent)).replace(os.sep, "/"), "relative_path": relative_path, "size_bytes": size, "sha256": digest})
    return outputs


def _all_outputs_payload(sample_root: Path) -> list[dict[str, object]]:
    files = [path for path in sorted(sample_root.rglob("*")) if path.is_file()]
    return _outputs_payload(
        sample_root,
        tuple(ExpectedOutput(str(path.relative_to(sample_root)).replace(os.sep, "/")) for path in files),
    )


def _verified_outputs_payload(
    sample_root: Path,
    expected: Sequence[ExpectedOutput],
    *,
    capture_extra_outputs: bool,
) -> list[dict[str, object]]:
    """Verify a declared output contract and return its published payload."""

    actual = {
        str(path.relative_to(sample_root)).replace(os.sep, "/")
        for path in sample_root.rglob("*")
        if path.is_file()
    } if sample_root.is_dir() else set()
    required = {item.relative_path.replace("\\", "/") for item in expected}
    missing = sorted(required - actual)
    if missing:
        raise OutputError(
            "expected batch output(s) are missing: " + ", ".join(missing)
        )
    unexpected = sorted(actual - required)
    if unexpected and not capture_extra_outputs:
        raise OutputError(
            "batch command produced undeclared output(s): " + ", ".join(unexpected)
        )
    selected = sorted(actual) if capture_extra_outputs else sorted(required)
    return _outputs_payload(
        sample_root,
        tuple(ExpectedOutput(relative_path) for relative_path in selected),
    )


def _resume_valid(
    job: dict[str, object],
    source: DiscoveredInput,
    spec: BatchCommandSpec,
    command: str,
    tail: Sequence[str],
    options: Namespace,
    run_dir: Path,
) -> bool:
    if job.get("status") not in {"succeeded", "threshold_failed", "failed"}:
        return False
    input_payload = job.get("input", {})
    if not isinstance(input_payload, dict):
        return False
    stored_identity = input_payload.get("source_identity") or _legacy_identity(input_payload.get("resolved_path"))
    if stored_identity not in {None, source.source_identity}:
        return False
    current_fingerprint = snapshot_file(source.resolved_path)[1]
    if input_payload.get("sha256") != current_fingerprint.sha256:
        return False
    outputs = job.get("outputs")
    if not isinstance(outputs, list):
        return False
    sample_root = run_dir / "outputs" / str(job.get("sample_key", ""))
    expected = spec.collect_outputs(sample_root, str(source.resolved_path), tail, options)
    stored_relative = {
        str(output.get("relative_path", "")).replace("\\", "/")
        for output in outputs
        if isinstance(output, dict)
    }
    expected_relative = {item.relative_path for item in expected}
    if job.get("status") == "failed":
        # Failed jobs deliberately publish no partial output.  Their verified
        # replay state is the source fingerprint plus the durable diagnostic
        # log, not a successful output contract.
        if outputs:
            return False
        if sample_root.exists() and any(sample_root.rglob("*")):
            return False
    else:
        if not expected_relative.issubset(stored_relative):
            return False
        if not outputs and spec.output_mode != "directory":
            return False
    actual_relative = {
        str(path.relative_to(sample_root)).replace(os.sep, "/")
        for path in sample_root.rglob("*")
        if path.is_file()
    }
    if job.get("status") != "failed" and stored_relative != actual_relative:
        return False
    if (
        job.get("status") != "failed"
        and not spec.capture_extra_outputs
        and actual_relative != expected_relative
    ):
        return False
    for output in outputs:
        if not isinstance(output, dict):
            return False
        if not _safe_relative(output.get("relative_path")) or not _safe_relative(output.get("path")):
            return False
        relative = Path(str(output["relative_path"]))
        path = run_dir / str(output.get("path", ""))
        if path != sample_root / relative:
            return False
        if not _within(path, sample_root):
            return False
        if not path.is_file():
            return False
        digest, size = _sha256(path)
        if digest != output.get("sha256") or size != output.get("size_bytes"):
            return False
    log = run_dir / str(job.get("stderr_log", ""))
    if not bool(job.get("stderr_log")) or not log.is_file():
        return False
    if not _within(log, run_dir):
        return False
    return hashlib.sha256(log.read_bytes()).hexdigest() == job.get("stderr_sha256")


def _run_one(
    spec: BatchCommandSpec,
    command: str,
    source: DiscoveredInput,
    sample_key: str,
    stage_root: Path,
    final_output_dir: Path,
    tail: Sequence[str],
    options: Namespace,
    input_payload: dict[str, object] | None = None,
) -> dict[str, object]:
    sample_root = stage_root / "outputs" / sample_key
    if spec.output_mode != "directory":
        sample_root.mkdir(parents=True, exist_ok=True)
    logs_root = stage_root / "logs"
    logs_root.mkdir(parents=True, exist_ok=True)
    expected = spec.collect_outputs(sample_root, str(source.resolved_path), tail, options)
    snapshot_root = Path(tempfile.mkdtemp(prefix=".gbparse-input.", dir=stage_root.parent))
    snapshot_path = snapshot_root / source.resolved_path.name
    started_clock = time.perf_counter()
    started = _now()
    fingerprint = None
    logical_argv = _logical_argv(
        spec,
        command,
        source,
        sample_key,
        tail,
        options,
        final_output_dir,
    )
    physical_argv = _build_argv(spec, command, snapshot_path, sample_root, tail, options)
    try:
        stderr = ""
        exit_code = 3
        try:
            raw, fingerprint = snapshot_file(source.resolved_path)
            snapshot_path.write_bytes(raw)
            physical_argv = _build_argv(
                spec, command, snapshot_path, sample_root, tail, options
            )
            child_environment = os.environ.copy()
            child_environment["GBPARSE_SOURCE_LABEL"] = source.display_path
            child_environment["GBPARSE_SOURCE_SHA256"] = fingerprint.sha256
            completed = subprocess.run(
                physical_argv,
                capture_output=True,
                text=True,
                shell=False,
                check=False,
                env=child_environment,
            )
            stderr = completed.stderr or ""
            if completed.stdout:
                stderr += "\n[child stdout captured]\n" + completed.stdout
            exit_code = completed.returncode
        except OSError as exc:
            stderr = str(exc)
            exit_code = 3

        log_name = f"{sample_key}.stderr.log"
        log_path = logs_root / log_name
        status = "succeeded" if exit_code == 0 else "threshold_failed" if exit_code == 1 else "failed"
        exit_class = {0: "success", 1: "validation_threshold", 2: "usage", 3: "input_or_job", 4: "output"}.get(exit_code, "unexpected")
        outputs: list[dict[str, object]] = []
        if exit_code in {0, 1}:
            try:
                outputs = _verified_outputs_payload(
                    sample_root,
                    expected,
                    capture_extra_outputs=spec.capture_extra_outputs,
                )
            except OutputError as exc:
                status = "failed"
                exit_code = 4
                exit_class = "output"
                stderr += f"\nERROR: {exc}\n"
        if status == "failed" and sample_root.exists():
            # A child may have created a partial file before failing. It is
            # never published as a plausible successful job artifact.
            shutil.rmtree(sample_root, ignore_errors=True)

        stderr = _scrub_diagnostics(
            stderr,
            stage_root=stage_root,
            snapshot_root=snapshot_root,
        )
        with log_path.open("w", encoding="utf-8", newline="") as log_handle:
            log_handle.write(stderr)
            log_handle.flush()
            os.fsync(log_handle.fileno())
        log_digest = hashlib.sha256(log_path.read_bytes()).hexdigest()
        result = {
            "input": _input_payload(source, fingerprint, fallback=input_payload),
            "status": status,
            "exit_code": exit_code,
            "exit_class": exit_class,
            # Persist the logical replay command, not the physical snapshot or
            # durable in-progress work-tree paths used by this execution.
            "argv": logical_argv,
            "started_at": started,
            "ended_at": _now(),
            "duration_seconds": round(time.perf_counter() - started_clock, 6),
            "outputs": outputs,
            "stderr_log": str(log_path.relative_to(stage_root)).replace(os.sep, "/"),
            "stderr_sha256": log_digest,
        }
        return result
    finally:
        # This must include subprocess interruption, output-contract failure,
        # and log-write failure.  The durable manifest is written by the
        # caller before entering this function and remains resumable when a
        # BaseException escapes.
        shutil.rmtree(snapshot_root, ignore_errors=True)


def _aggregate(manifest: dict[str, object]) -> tuple[int, int]:
    jobs = manifest.get("jobs", [])
    failed = sum(
        1
        for job in jobs
        if job.get("status") in {"failed", "threshold_failed"}
    )
    codes = [
        int(job.get("exit_code", 0) or 0)
        for job in jobs
        if job.get("status") in {"failed", "threshold_failed", "skipped_unchanged"}
    ]
    if any(code == 4 for code in codes):
        return 4, failed
    if any(code == 3 for code in codes):
        return 3, failed
    if any(code == 2 for code in codes):
        return 2, failed
    if any(code == 1 for code in codes):
        return 1, failed
    if any(code not in {0, 1, 2, 3, 4} or code < 0 for code in codes):
        return 3, failed
    return 0, failed


def _legacy_identity(value: object) -> str | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        return os.path.normcase(os.fspath(Path(value).resolve(strict=False)))
    except OSError:
        return os.path.normcase(os.path.abspath(value))


def _safe_relative(value: object) -> bool:
    if not isinstance(value, str) or not value:
        return False
    path = Path(value)
    return not path.is_absolute() and ".." not in path.parts


def _safe_sample_key(value: object) -> bool:
    if not isinstance(value, str) or not value:
        return False
    path = Path(value)
    return not path.is_absolute() and len(path.parts) == 1 and path.parts[0] not in {".", ".."}


def _within(path: Path, root: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(root.resolve(strict=False))
    except ValueError:
        return False
    return True


def _validate_manifest(manifest: object) -> dict[str, object]:
    """Validate the persisted v1 manifest before trusting any resume state."""

    try:
        import jsonschema

        schema = json.loads(
            resources.files("genbank_parser")
            .joinpath("data", "schemas", "batch-v1.schema.json")
            .read_text(encoding="utf-8")
        )
        jsonschema.validate(manifest, schema)
    except ModuleNotFoundError:
        pass
    except (OSError, json.JSONDecodeError, jsonschema.ValidationError) as exc:  # type: ignore[union-attr]
        raise InputError(f"batch manifest does not validate against batch-v1 schema: {exc}") from exc
    if not isinstance(manifest, dict):
        raise InputError("batch manifest must be a JSON object")
    required = {
        "schema_version",
        "gbparse_version",
        "python_version",
        "biopython_version",
        "platform",
        "created_at",
        "updated_at",
        "command",
        "command_args",
        "runner",
        "jobs",
    }
    if not required.issubset(manifest):
        missing = sorted(required - set(manifest))
        raise InputError(f"batch manifest is missing required fields: {', '.join(missing)}")
    if manifest.get("schema_version") != "gbparse.batch.v1":
        raise InputError("incompatible batch manifest schema")
    if not isinstance(manifest.get("command_args"), list) or not all(
        isinstance(value, str) for value in manifest["command_args"]
    ):
        raise InputError("batch manifest command_args must be a string array")
    runner = manifest.get("runner")
    if not isinstance(runner, dict) or not isinstance(runner.get("jobs"), int) or runner["jobs"] <= 0:
        raise InputError("batch manifest runner is invalid")
    working_directory = runner.get("working_directory")
    if working_directory is None:
        if manifest.get("gbparse_version") == __version__:
            raise InputError(
                "new batch manifests must record runner.working_directory"
            )
    elif not isinstance(working_directory, str) or not os.path.isabs(working_directory):
        raise InputError("batch manifest runner.working_directory must be absolute")
    if not isinstance(manifest.get("jobs"), list):
        raise InputError("batch manifest jobs must be an array")
    statuses = {"pending", "running", "succeeded", "threshold_failed", "failed", "skipped_unchanged", "not_requested"}
    job_ids: set[str] = set()
    sample_keys: set[str] = set()
    source_identities: set[str] = set()
    output_paths: set[str] = set()
    for job in manifest["jobs"]:
        if not isinstance(job, dict):
            raise InputError("batch manifest contains a non-object job")
        for key in ("job_id", "sample_key", "status", "input"):
            if key not in job:
                raise InputError(f"batch manifest job is missing {key!r}")
        job_id = str(job["job_id"])
        if job_id in job_ids:
            raise InputError(f"batch manifest has duplicate job_id {job_id!r}")
        job_ids.add(job_id)
        status = str(job["status"])
        if status not in statuses:
            raise InputError(f"batch manifest has invalid job status {status!r}")
        if "resume_action" in job and job["resume_action"] not in {
            "reused_unchanged",
            "rerun_changed",
            "replayed_running",
            "deferred_after_stop",
            "not_requested",
        }:
            raise InputError("batch manifest has an invalid resume_action")
        input_payload = job["input"]
        if not isinstance(input_payload, dict):
            raise InputError("batch manifest job input must be an object")
        if status != "not_requested":
            source_identity = input_payload.get("source_identity") or _legacy_identity(
                input_payload.get("resolved_path")
            )
            if source_identity is not None:
                source_identity = str(source_identity)
                if source_identity in source_identities:
                    raise InputError(
                        f"batch manifest has duplicate active source identity {source_identity!r}"
                    )
                source_identities.add(source_identity)
            sample_key = str(job["sample_key"])
            if not _safe_sample_key(sample_key):
                raise InputError(f"batch manifest has unsafe sample_key {sample_key!r}")
            sample_key_identity = os.path.normcase(sample_key)
            if sample_key_identity in sample_keys:
                raise InputError(f"batch manifest has duplicate active sample_key {sample_key!r}")
            sample_keys.add(sample_key_identity)
        outputs = job.get("outputs", [])
        if not isinstance(outputs, list):
            raise InputError("batch manifest job outputs must be an array")
        for output in outputs:
            if not isinstance(output, dict):
                raise InputError("batch manifest output must be an object")
            if not {"path", "relative_path", "size_bytes", "sha256"}.issubset(output):
                raise InputError("batch manifest output is missing required fields")
            if not _safe_relative(output.get("relative_path")) or not _safe_relative(output.get("path")):
                raise InputError("batch manifest contains an unsafe output path")
            path = str(output.get("path", ""))
            if status != "not_requested":
                path_identity = os.path.normcase(path.replace("\\", "/"))
                if path_identity in output_paths:
                    raise InputError(f"batch manifest has duplicate active output path {path!r}")
                output_paths.add(path_identity)
        if "stderr_log" in job and not _safe_relative(job.get("stderr_log")):
            raise InputError("batch manifest contains an unsafe stderr log path")
    return manifest


def _normalize_legacy_outcomes(manifest: dict[str, object]) -> None:
    """Load pre-hardening manifests without retaining a fake terminal status."""

    for job in manifest.get("jobs", []):
        if not isinstance(job, dict) or job.get("status") != "skipped_unchanged":
            continue
        code = job.get("exit_code")
        if code == 1:
            job["status"] = "threshold_failed"
        elif code == 0:
            job["status"] = "succeeded"
        else:
            job["status"] = "failed"
        job["resume_action"] = "reused_unchanged"


def _remove_job_artifacts(run_dir: Path, job: dict[str, object], sample_key: str) -> None:
    """Remove only one validated job's output and log artifacts."""

    if not _safe_sample_key(sample_key):
        return
    output_root = run_dir / "outputs" / sample_key
    if output_root.is_dir():
        shutil.rmtree(output_root)
    stderr_log = job.get("stderr_log")
    if _safe_relative(stderr_log):
        log_path = run_dir / str(stderr_log)
        if _within(log_path, run_dir) and log_path.is_file():
            log_path.unlink()


def _reset_job_for_replay(
    job: dict[str, object],
    source: DiscoveredInput,
    sample_key: str,
    *,
    action: str | None,
) -> None:
    """Keep identity while clearing terminal execution state."""

    job.update(
        {
            "job_id": _job_id(source),
            "sample_key": sample_key,
            "status": "pending",
            "input": _input_payload(source),
            "outputs": [],
        }
    )
    if action is None:
        job.pop("resume_action", None)
    else:
        job["resume_action"] = action
    for key in (
        "exit_code",
        "exit_class",
        "argv",
        "stderr_log",
        "stderr_sha256",
        "started_at",
        "ended_at",
        "duration_seconds",
    ):
        job.pop(key, None)


def _inprogress_path(output_path: Path) -> Path:
    return output_path.parent / f".{output_path.name}.gbparse-inprogress"


def _remove_path(path: Path) -> None:
    """Remove either a directory tree or a file/symlink safely."""

    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def execute_batch(
    inputs: Iterable[str | Path],
    *,
    command: str,
    output_dir: str | Path,
    tail: Sequence[str] = (),
    jobs: int = 1,
    on_error: str = "continue",
    resume: bool = False,
    force: bool = False,
) -> BatchExecutionResult:
    """Execute a registered command with atomic per-job publication."""

    inputs = tuple(inputs)
    if tail and tail[0] == "--":
        tail = tuple(tail[1:])
    if jobs <= 0:
        raise ValueError("jobs must be positive")
    if on_error not in {"continue", "stop"}:
        raise ValueError("on_error must be continue or stop")
    if resume and force:
        raise ValueError("--resume and --force are mutually exclusive")
    spec = get_command_spec(command, tail)
    _validate_tail(spec, tail)
    _validate_batch_dependencies(command, tail)
    sentinel = Path("__gbparse_batch_input__.gbff")
    base_options = _validate_command_argv(command, [str(sentinel), *tail])
    raw_validation = spec.build_job_argv(command, str(sentinel), tail, base_options)
    validation_argv = [value.replace("{OUTPUT}", "__gbparse_batch_output__") for value in raw_validation[1:]]
    options = _validate_command_argv(command, validation_argv)
    output_path = Path(output_dir)
    reject_input_output_collision(inputs, output_path)
    inprogress_path = _inprogress_path(output_path)
    staging_siblings = (
        tuple(output_path.parent.glob(f".{output_path.name}.*"))
        if output_path.parent.exists()
        else ()
    )
    discovered = discover_inputs(
        inputs,
        exclude=(output_path, inprogress_path, *staging_siblings),
    )
    if not discovered:
        raise InputError("no GenBank inputs were discovered")
    keys = assign_sample_keys(discovered)

    if resume:
        if inprogress_path.is_dir() and _manifest_path(inprogress_path).is_file():
            resume_dir = inprogress_path
        elif output_path.is_dir() and _manifest_path(output_path).is_file():
            resume_dir = output_path
        else:
            raise InputError(
                "--resume requires a published run or durable in-progress batch manifest"
            )
        try:
            manifest = json.loads(
                _manifest_path(resume_dir).read_text(encoding="utf-8")
            )
        except (OSError, json.JSONDecodeError) as exc:
            raise InputError(f"could not read batch manifest: {exc}") from exc
        manifest = _validate_manifest(manifest)
        _normalize_legacy_outcomes(manifest)
        if (
            manifest.get("command") != command
            or manifest.get("command_args") != list(tail)
            or manifest.get("gbparse_version") != __version__
            or manifest.get("python_version") != platform.python_version()
            or manifest.get("biopython_version") != __import__("Bio").__version__
            or manifest.get("platform") != platform.platform()
        ):
            raise InputError("batch resume manifest environment or command does not match")
        runner = manifest["runner"]
        if runner.get("jobs") != jobs or runner.get("on_error") != on_error:
            raise InputError("batch resume runner settings do not match the original run")
        if resume_dir == output_path:
            if inprogress_path.exists():
                raise InputError(
                    f"a durable in-progress batch already exists; resume it directly: {inprogress_path}"
                )
            resume_copy = output_path.parent / (
                f".{output_path.name}.gbparse-resume-copy.{secrets.token_hex(8)}"
            )
            try:
                shutil.copytree(output_path, resume_copy)
                os.replace(resume_copy, inprogress_path)
            except BaseException as exc:
                cleanup_errors: list[OSError] = []
                for candidate in (resume_copy, inprogress_path):
                    if not (candidate.exists() or candidate.is_symlink()):
                        continue
                    try:
                        _remove_path(candidate)
                    except OSError as cleanup_exc:
                        cleanup_errors.append(cleanup_exc)
                if cleanup_errors:
                    retained = tuple(
                        candidate
                        for candidate in (resume_copy, inprogress_path)
                        if candidate.exists() or candidate.is_symlink()
                    )
                    raise OutputRecoveryError(
                        "could not clean partial batch resume copy; published run remains authoritative",
                        destination=output_path,
                        staging=retained,
                    ) from cleanup_errors[0]
                if isinstance(exc, OSError):
                    raise OutputError(
                        f"could not create durable batch resume tree {inprogress_path}: {exc}"
                    ) from exc
                raise
        work_dir = inprogress_path
    else:
        if output_path.exists() and not force:
            raise OutputError(f"batch output directory already exists; use --force: {output_path}")
        if inprogress_path.exists():
            raise InputError(
                f"durable in-progress batch exists; use --resume or remove it after review: {inprogress_path}"
            )
        manifest = _manifest_template(command, tail, discovered, keys, jobs, on_error)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            inprogress_path.mkdir()
            (inprogress_path / "outputs").mkdir()
            (inprogress_path / "logs").mkdir()
        except OSError as exc:
            raise OutputError(
                f"could not create durable batch work tree {inprogress_path}: {exc}"
            ) from exc
        work_dir = inprogress_path
        _write_manifest(_manifest_path(work_dir), manifest)

    for job in manifest.get("jobs", []):
        if job.get("status") == "running":
            _remove_job_artifacts(work_dir, job, str(job.get("sample_key", "")))
            job["status"] = "pending"
            job["resume_action"] = "replayed_running"
            job["outputs"] = []
            for key in (
                "exit_code",
                "exit_class",
                "argv",
                "stderr_log",
                "stderr_sha256",
                "started_at",
                "ended_at",
                "duration_seconds",
            ):
                job.pop(key, None)

    existing_by_identity: dict[str, dict[str, object]] = {}
    for job in manifest.get("jobs", []):
        payload = job.get("input", {})
        identity = payload.get("source_identity") if isinstance(payload, dict) else None
        identity = identity or _legacy_identity(
            payload.get("resolved_path") if isinstance(payload, dict) else None
        )
        if identity is not None:
            existing_by_identity[str(identity)] = job

    current_identities = {source.source_identity for source in discovered}
    for job in manifest.get("jobs", []):
        payload = job.get("input", {})
        identity = payload.get("source_identity") if isinstance(payload, dict) else None
        identity = identity or _legacy_identity(
            payload.get("resolved_path") if isinstance(payload, dict) else None
        )
        if identity not in current_identities:
            old_key = str(job.get("sample_key", ""))
            if job.get("status") != "not_requested":
                _remove_job_artifacts(work_dir, job, old_key)
            job["status"] = "not_requested"
            job["resume_action"] = "not_requested"
            job["outputs"] = []
            for key in (
                "exit_code",
                "exit_class",
                "argv",
                "stderr_log",
                "stderr_sha256",
                "started_at",
                "ended_at",
                "duration_seconds",
            ):
                job.pop(key, None)

    pending: list[tuple[DiscoveredInput, str, dict[str, object]]] = []
    stop_barrier_position: int | None = None
    for position, source in enumerate(discovered):
        key = keys[source]
        existing = existing_by_identity.get(source.source_identity)
        if (
            resume
            and existing is not None
            and str(existing.get("sample_key", "")) == key
            and _resume_valid(existing, source, spec, command, tail, options, work_dir)
        ):
            if existing.get("status") == "skipped_unchanged":
                existing["status"] = (
                    "threshold_failed"
                    if existing.get("exit_code") == 1
                    else "succeeded"
                )
            existing["resume_action"] = "reused_unchanged"
            if on_error == "stop" and existing.get("status") in {
                "failed",
                "threshold_failed",
            }:
                if stop_barrier_position is None:
                    stop_barrier_position = position
            input_payload = existing.get("input")
            if isinstance(input_payload, dict):
                input_payload.update(
                    {
                        "requested_path": source.requested_path,
                        "source_identity": source.source_identity,
                        "display_path": source.display_path,
                        "resolved_path": str(source.resolved_path),
                    }
                )
            continue
        if existing is None:
            existing = {
                "job_id": _job_id(source),
                "sample_key": key,
                "status": "pending",
                "input": _input_payload(source),
                "outputs": [],
            }
            manifest.setdefault("jobs", []).append(existing)
            existing_by_identity[source.source_identity] = existing
        else:
            old_sample = str(existing.get("sample_key", key))
            old_status = str(existing.get("status", "pending"))
            prior_resume_action = existing.get("resume_action")
            _remove_job_artifacts(work_dir, existing, old_sample)
            replay_action = (
                "replayed_running"
                if old_status == "running" or prior_resume_action == "replayed_running"
                else "rerun_changed"
                if old_status in {"succeeded", "threshold_failed", "failed"}
                else None
            )
            _reset_job_for_replay(
                existing,
                source,
                key,
                action=replay_action,
            )
        sample_root = work_dir / "outputs" / key
        if sample_root.exists():
            shutil.rmtree(sample_root)
        if (
            on_error == "stop"
            and stop_barrier_position is not None
            and position > stop_barrier_position
        ):
            existing["resume_action"] = "deferred_after_stop"
            continue
        pending.append((source, key, existing))

    _validate_manifest(manifest)
    _write_manifest(_manifest_path(work_dir), manifest)

    def record_result(
        source: DiscoveredInput,
        existing: dict[str, object],
        result: dict[str, object],
    ) -> bool:
        resume_action = existing.get("resume_action")
        result_input = result.get("input")
        if not isinstance(result_input, dict):
            prior_input = existing.get("input")
            result_input = _input_payload(
                source,
                fallback=prior_input if isinstance(prior_input, dict) else None,
            )
        existing.update({**result, "input": result_input})
        if resume_action is not None:
            existing["resume_action"] = resume_action
        manifest["updated_at"] = _now()
        _write_manifest(_manifest_path(work_dir), manifest)
        return result["status"] in {"failed", "threshold_failed"}

    if pending:
        if jobs > 1 and on_error == "continue" and len(pending) > 1:
            for _source, _key, existing in pending:
                existing["status"] = "running"
            _write_manifest(_manifest_path(work_dir), manifest)
            with ThreadPoolExecutor(max_workers=min(jobs, len(pending))) as executor:
                futures = {
                    executor.submit(
                        _run_one,
                        spec,
                        command,
                        source,
                        key,
                        work_dir,
                        output_path,
                        tail,
                        options,
                        existing.get("input") if isinstance(existing.get("input"), dict) else None,
                    ): (source, existing)
                    for source, key, existing in pending
                }
                for future in as_completed(futures):
                    source, existing = futures[future]
                    record_result(source, existing, future.result())
        else:
            for source, key, existing in pending:
                existing["status"] = "running"
                _write_manifest(_manifest_path(work_dir), manifest)
                result = _run_one(
                    spec,
                    command,
                    source,
                    key,
                    work_dir,
                    output_path,
                    tail,
                    options,
                    existing.get("input") if isinstance(existing.get("input"), dict) else None,
                )
                if record_result(source, existing, result) and on_error == "stop":
                    stop_latched = True
                    break

    manifest["updated_at"] = _now()
    _validate_manifest(manifest)
    _write_manifest(_manifest_path(work_dir), manifest)
    # A successful directory move is the publication commit point.  Any
    # exception before it deliberately leaves the durable in-progress tree.
    # Auxiliary resources are rejected during preflight in v0.9.3.  Keep the
    # tuple boundary explicit so a future dependency snapshot can share the
    # same collision and publication guard without weakening this caller.
    auxiliary_inputs: tuple[str | Path, ...] = ()
    publish_directory_tree(
        work_dir,
        output_path,
        force=force or resume,
        inputs=(*inputs, *auxiliary_inputs),
        retain_staging_on_failure=True,
    )
    exit_code, failed_count = _aggregate(manifest)
    return BatchExecutionResult(manifest=manifest, failed_count=failed_count, exit_code=exit_code)


__all__ = [
    "BatchCommandSpec",
    "BatchExecutionResult",
    "BatchUsageError",
    "ExpectedOutput",
    "execute_batch",
    "get_command_spec",
    "registered_commands",
]
