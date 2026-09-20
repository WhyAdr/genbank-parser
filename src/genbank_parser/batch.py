"""Manifest-driven execution of registered single-input gbparse commands."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

from . import __version__
from .cli_io import (
    InputError,
    OutputError,
    publish_directory_tree,
    reject_input_output_collision,
    write_text,
)
from .discovery import (
    DiscoveredInput,
    assign_sample_keys,
    discover_inputs,
    fingerprint_file,
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
    build_job_argv: Callable[[str, str, Sequence[str]], list[str]] = field(
        default=lambda command, input_path, tail: [command, input_path, *tail]
    )
    collect_outputs: Callable[[Path, str, Sequence[str]], tuple[ExpectedOutput, ...]] = field(
        default=lambda sample_dir, input_path, tail: (ExpectedOutput("result.out"),)
    )


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


def _format_suffix(tail: Sequence[str], fallback: str) -> str:
    try:
        index = next(index for index, value in enumerate(tail) if value == "--format")
        value = tail[index + 1]
    except (StopIteration, IndexError):
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


def _source_stem(input_path: str) -> str:
    name = Path(input_path).name
    if name.casefold().endswith(".gz"):
        name = Path(name).stem
    return Path(name).stem


def _file_spec(name: str, suffix: str) -> BatchCommandSpec:
    def build(command: str, input_path: str, tail: Sequence[str]) -> list[str]:
        output_name = f"result{_format_suffix(tail, suffix)}"
        return [command, input_path, *tail, "--output", "{OUTPUT}/" + output_name]

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str]) -> tuple[ExpectedOutput, ...]:
        return (ExpectedOutput(f"result{_format_suffix(tail, suffix)}"),)

    return BatchCommandSpec(name=name, default_suffix=suffix, build_job_argv=build, collect_outputs=collect)


def _positional_spec(name: str, suffix: str) -> BatchCommandSpec:
    def build(command: str, input_path: str, tail: Sequence[str]) -> list[str]:
        return [command, input_path, *tail, "{OUTPUT}/result" + _format_suffix(tail, suffix)]

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str]) -> tuple[ExpectedOutput, ...]:
        return (ExpectedOutput("result" + _format_suffix(tail, suffix)),)

    return BatchCommandSpec(
        name=name,
        default_suffix=suffix,
        build_job_argv=build,
        collect_outputs=collect,
        forbidden_options=("--output", "--output-dir", "--force"),
    )


def _sequence_spec() -> BatchCommandSpec:
    def build(command: str, input_path: str, tail: Sequence[str]) -> list[str]:
        return [command, input_path, *tail, "--fna", "{OUTPUT}/genome.fna", "--ffn", "{OUTPUT}/cds.ffn"]

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str]) -> tuple[ExpectedOutput, ...]:
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
    def build(command: str, input_path: str, tail: Sequence[str]) -> list[str]:
        argv = [command, input_path, *tail, "--output", "{OUTPUT}/result" + _format_suffix(tail, ".json")]
        if "--visualize" in tail or "--viz-output" in tail:
            argv.extend(("--visualize", "--viz-output", "{OUTPUT}/neighborhood.svg"))
        return argv

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str]) -> tuple[ExpectedOutput, ...]:
        result = [ExpectedOutput("result" + _format_suffix(tail, ".json"))]
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
    def build(command: str, input_path: str, tail: Sequence[str]) -> list[str]:
        return [command, input_path, *tail, "--output", "{OUTPUT}/report" + _format_suffix(tail, ".json"), "--output-dir", "{OUTPUT}/markers"]

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str]) -> tuple[ExpectedOutput, ...]:
        return (ExpectedOutput("report" + _format_suffix(tail, ".json")),)

    return BatchCommandSpec(
        name="phylo",
        output_mode="multi_file",
        default_suffix=None,
        forbidden_options=("--output", "--output-dir", "--force"),
        build_job_argv=build,
        collect_outputs=collect,
    )


def _export_spec() -> BatchCommandSpec:
    def build(command: str, input_path: str, tail: Sequence[str]) -> list[str]:
        if "ncbi-table" in tail:
            return [command, input_path, *tail, "--output-dir", "{OUTPUT}"]
        return [command, input_path, *tail, "--output", "{OUTPUT}/result" + _format_suffix(tail, ".out")]

    def collect(sample_dir: Path, input_path: str, tail: Sequence[str]) -> tuple[ExpectedOutput, ...]:
        if "ncbi-table" in tail:
            return tuple(ExpectedOutput(str(path.relative_to(sample_dir))) for path in sorted(sample_dir.rglob("*")) if path.is_file())
        return (ExpectedOutput("result" + _format_suffix(tail, ".out")),)

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
    if command == "export" and "ncbi-table" in tail:
        return replace(spec, output_mode="directory")
    return spec


def _validate_tail(spec: BatchCommandSpec, tail: Sequence[str]) -> None:
    for token in tail:
        for forbidden in spec.forbidden_options:
            if token == forbidden or token.startswith(forbidden + "="):
                raise BatchUsageError(
                    f"batch owns {forbidden}; remove it from the wrapped command tail"
                )


def _validate_command_argv(command: str, argv: Sequence[str]) -> None:
    from .cli import create_command_parser

    parser = create_command_parser(command)
    try:
        parser.parse_args(list(argv))
    except SystemExit as exc:
        raise BatchUsageError(
            f"batch command {command!r} arguments failed pre-validation"
        ) from exc


def _build_argv(spec: BatchCommandSpec, command: str, input_path: Path, sample_dir: Path, tail: Sequence[str]) -> list[str]:
    raw = spec.build_job_argv(command, str(input_path), tail)
    return [
        sys.executable,
        "-m",
        "genbank_parser.cli",
        *[str(sample_dir) if value == "{OUTPUT}" else str(sample_dir / value[len("{OUTPUT}/") :]) if value.startswith("{OUTPUT}/") else value for value in raw],
    ]


def _manifest_path(run_dir: Path) -> Path:
    return run_dir / "batch-manifest.json"


def _write_manifest(path: Path, manifest: dict[str, object]) -> None:
    write_text(json.dumps(manifest, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2) + "\n", path, force=True)


def _job_id(source: DiscoveredInput) -> str:
    return "job-" + hashlib.sha256(source.display_path.encode("utf-8")).hexdigest()[:16]


def _input_payload(source: DiscoveredInput) -> dict[str, object]:
    fp = fingerprint_file(source.resolved_path)
    return {
        "requested_path": source.requested_path,
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
        "runner": {"jobs": jobs, "on_error": on_error, "discovery_roots": sorted({str(item.discovery_root) for item in inputs if item.discovery_root is not None})},
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
        outputs.append({"path": str(path.relative_to(sample_root.parent.parent)).replace(os.sep, "/"), "relative_path": item.relative_path, "size_bytes": size, "sha256": digest})
    return outputs


def _all_outputs_payload(sample_root: Path) -> list[dict[str, object]]:
    files = [path for path in sorted(sample_root.rglob("*")) if path.is_file()]
    return _outputs_payload(sample_root, tuple(ExpectedOutput(str(path.relative_to(sample_root))) for path in files))


def _resume_valid(job: dict[str, object], source: DiscoveredInput, command: str, tail: Sequence[str], run_dir: Path) -> bool:
    if job.get("status") not in {"succeeded", "threshold_failed"}:
        return False
    if job.get("input", {}).get("display_path") != source.display_path:
        return False
    if job.get("input", {}).get("sha256") != fingerprint_file(source.resolved_path).sha256:
        return False
    for output in job.get("outputs", []):
        path = run_dir / str(output.get("path", ""))
        if not path.is_file():
            return False
        digest, size = _sha256(path)
        if digest != output.get("sha256") or size != output.get("size_bytes"):
            return False
    log = run_dir / str(job.get("stderr_log", ""))
    return bool(job.get("stderr_log")) and log.is_file()


def _run_one(
    spec: BatchCommandSpec,
    command: str,
    source: DiscoveredInput,
    sample_key: str,
    stage_root: Path,
    tail: Sequence[str],
) -> dict[str, object]:
    sample_root = stage_root / "outputs" / sample_key
    if not (command == "export" and "ncbi-table" in tail):
        sample_root.mkdir(parents=True, exist_ok=True)
    logs_root = stage_root / "logs"
    logs_root.mkdir(parents=True, exist_ok=True)
    expected = spec.collect_outputs(sample_root, str(source.resolved_path), tail)
    argv = _build_argv(spec, command, source.resolved_path, sample_root, tail)
    started_clock = time.perf_counter()
    started = _now()
    try:
        completed = subprocess.run(argv, capture_output=True, text=True, shell=False, check=False)
        stderr = completed.stderr
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
            outputs = (
                _all_outputs_payload(sample_root)
                if (command == "export" and "ncbi-table" in tail) or command == "phylo"
                else _outputs_payload(sample_root, expected)
            )
        except OutputError as exc:
            status = "failed"
            exit_code = 4
            exit_class = "output"
            stderr += f"\nERROR: {exc}\n"
    if status == "failed" and sample_root.exists():
        # A child may have created a partial file before failing. It is never
        # published as a plausible successful job artifact.
        shutil.rmtree(sample_root, ignore_errors=True)
    log_path.write_text(stderr, encoding="utf-8")
    return {
        "status": status,
        "exit_code": exit_code,
        "exit_class": exit_class,
        "argv": argv,
        "started_at": started,
        "ended_at": _now(),
        "duration_seconds": round(time.perf_counter() - started_clock, 6),
        "outputs": outputs,
        "stderr_log": str(log_path.relative_to(stage_root)).replace(os.sep, "/"),
        "stderr_sha256": hashlib.sha256(stderr.encode("utf-8")).hexdigest(),
    }


def _aggregate(manifest: dict[str, object]) -> tuple[int, int]:
    jobs = manifest.get("jobs", [])
    failed = sum(1 for job in jobs if job.get("status") == "failed")
    codes = [int(job.get("exit_code", 0) or 0) for job in jobs if job.get("status") in {"failed", "threshold_failed"}]
    if any(code == 4 for code in codes):
        return 4, failed
    if any(code == 3 for code in codes):
        return 3, failed
    if any(code == 2 for code in codes):
        return 2, failed
    if any(code == 1 for code in codes):
        return 1, failed
    return 0, failed


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
    output_path = Path(output_dir)
    reject_input_output_collision(inputs, output_path)
    staging_siblings = tuple(output_path.parent.glob(f".{output_path.name}.*")) if output_path.parent.exists() else ()
    discovered = discover_inputs(inputs, exclude=(output_path, *staging_siblings))
    if not discovered:
        raise InputError("no GenBank inputs were discovered")
    keys = assign_sample_keys(discovered)
    # Parse once before any child starts. The injected paths are parser-valid
    # sentinels; biological execution remains isolated in child processes.
    sentinel = Path("__gbparse_batch_input__.gbff")
    sample_sentinel = Path("__gbparse_batch_output__")
    validation_argv = [str(sentinel)]
    raw_validation = spec.build_job_argv(command, str(sentinel), tail)
    validation_argv = [value.replace("{OUTPUT}", str(sample_sentinel)) for value in raw_validation[1:]]
    _validate_command_argv(command, validation_argv)

    if resume:
        manifest_file = _manifest_path(output_path)
        if not output_path.is_dir() or not manifest_file.is_file():
            raise InputError("--resume requires an existing batch run directory and manifest")
        try:
            manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise InputError(f"could not read batch manifest: {exc}") from exc
        if manifest.get("schema_version") != "gbparse.batch.v1":
            raise InputError("incompatible batch manifest schema")
        if manifest.get("command") != command or manifest.get("command_args") != list(tail) or manifest.get("gbparse_version") != __version__:
            raise InputError("batch resume manifest command or gbparse version does not match")
    else:
        if output_path.exists() and not force:
            raise OutputError(f"batch output directory already exists; use --force: {output_path}")
        manifest = _manifest_template(command, tail, discovered, keys, jobs, on_error)

    parent = output_path.parent
    parent.mkdir(parents=True, exist_ok=True)
    work_dir = Path(tempfile.mkdtemp(prefix=f".{output_path.name}.work.", dir=parent))
    try:
        if resume:
            shutil.copytree(output_path, work_dir, dirs_exist_ok=True)
            for job in manifest.get("jobs", []):
                if job.get("status") == "running":
                    job["status"] = "pending"
        else:
            (work_dir / "outputs").mkdir()
            (work_dir / "logs").mkdir()
        existing_by_display = {
            str(job.get("input", {}).get("display_path")): job for job in manifest.get("jobs", [])
        }
        current_display = {source.display_path for source in discovered}
        for job in manifest.get("jobs", []):
            if job.get("input", {}).get("display_path") not in current_display:
                job["status"] = "not_requested"
        pending: list[tuple[DiscoveredInput, str, dict[str, object]]] = []
        for source in discovered:
            key = keys[source]
            existing = existing_by_display.get(source.display_path)
            if resume and existing is not None and _resume_valid(existing, source, command, tail, output_path):
                existing["status"] = "skipped_unchanged"
                continue
            if existing is None:
                existing = {"job_id": _job_id(source), "sample_key": key}
                manifest.setdefault("jobs", []).append(existing)
            else:
                old_sample = str(existing.get("sample_key", key))
                if old_sample != key:
                    existing["sample_key"] = key
            sample_root = work_dir / "outputs" / key
            if sample_root.exists():
                shutil.rmtree(sample_root)
            pending.append((source, key, existing))

        def record_result(
            source: DiscoveredInput,
            existing: dict[str, object],
            result: dict[str, object],
        ) -> bool:
            existing.update({"input": _input_payload(source), **result})
            manifest["updated_at"] = _now()
            _write_manifest(_manifest_path(work_dir), manifest)
            return result["status"] == "failed"

        if jobs > 1 and on_error == "continue" and len(pending) > 1:
            with ThreadPoolExecutor(max_workers=min(jobs, len(pending))) as executor:
                futures = {
                    executor.submit(_run_one, spec, command, source, key, work_dir, tail): (source, existing)
                    for source, key, existing in pending
                }
                for future in as_completed(futures):
                    source, existing = futures[future]
                    record_result(source, existing, future.result())
        else:
            for source, key, existing in pending:
                result = _run_one(spec, command, source, key, work_dir, tail)
                if record_result(source, existing, result) and on_error == "stop":
                    break
        manifest["updated_at"] = _now()
        _write_manifest(_manifest_path(work_dir), manifest)
        # A final manifest exists before the directory move, so an interrupted
        # publication never presents a successful artifact without provenance.
        publish_directory_tree(work_dir, output_path, force=force or resume, inputs=inputs)
        work_dir = None  # type: ignore[assignment]
    finally:
        if work_dir is not None and work_dir.exists():
            shutil.rmtree(work_dir, ignore_errors=True)
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
