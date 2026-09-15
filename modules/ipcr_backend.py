"""Resolve and build the ipcr backend used by probe-tester.

The preferred source is the pinned Git submodule at ``third_party/ipcr``.
When that source is available, probe-tester builds a native executable into
``bin/`` and records the source revision beside it. An explicit ``--ipcr-bin``
override or a system ``ipcr`` executable can still be used.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Optional

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
IPCR_SOURCE_DIR = PROJECT_ROOT / "third_party" / "ipcr"
IPCR_BIN_DIR = PROJECT_ROOT / "bin"
IPCR_BINARY = IPCR_BIN_DIR / ("ipcr.exe" if sys.platform.startswith("win") else "ipcr")
IPCR_REVISION_FILE = IPCR_BIN_DIR / "ipcr.source-revision"


class IPCRBackendError(RuntimeError):
    """Raised when the ipcr executable cannot be resolved or built."""


def _source_is_available(source_dir: Path = IPCR_SOURCE_DIR) -> bool:
    return (source_dir / "go.mod").is_file() and (
        source_dir / "cmd" / "ipcr" / "main.go"
    ).is_file()


def _git_text(source_dir: Path, *args: str) -> Optional[str]:
    git = shutil.which("git")
    if not git:
        return None
    try:
        result = subprocess.run(
            [git, "-C", str(source_dir), *args],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    value = result.stdout.strip()
    return value or None


def _source_revision(source_dir: Path = IPCR_SOURCE_DIR) -> Optional[str]:
    return _git_text(source_dir, "rev-parse", "HEAD")


def _source_version(source_dir: Path = IPCR_SOURCE_DIR) -> str:
    return (
        _git_text(source_dir, "describe", "--tags", "--dirty", "--always")
        or "submodule"
    )


def _binary_matches_revision(binary: Path, revision: Optional[str]) -> bool:
    if not binary.is_file():
        return False
    if revision is None:
        # Source archives without .git metadata can still use an existing build.
        return True
    try:
        return IPCR_REVISION_FILE.read_text(encoding="utf-8").strip() == revision
    except OSError:
        return False


def build_ipcr_from_submodule(*, force: bool = False) -> Path:
    """Build the pinned ``third_party/ipcr`` source and return its executable."""
    if not _source_is_available():
        raise IPCRBackendError(
            "ipcr submodule is not initialized at third_party/ipcr. "
            "Run 'git submodule update --init --recursive'."
        )

    revision = _source_revision()
    if not force and _binary_matches_revision(IPCR_BINARY, revision):
        return IPCR_BINARY

    go = shutil.which("go")
    if not go:
        raise IPCRBackendError(
            "Go is required to build the ipcr submodule (ipcr requires Go >= 1.22), "
            "or provide --ipcr-bin /path/to/ipcr."
        )

    IPCR_BIN_DIR.mkdir(parents=True, exist_ok=True)
    tmp_binary = IPCR_BINARY.with_name(IPCR_BINARY.name + ".tmp")
    version = _source_version()
    cmd = [
        go,
        "build",
        "-trimpath",
        "-ldflags",
        f"-s -w -X ipcr/internal/version.Version={version}",
        "-o",
        str(tmp_binary),
        "./cmd/ipcr",
    ]

    logger.info("Building ipcr %s from %s", version, IPCR_SOURCE_DIR)
    try:
        result = subprocess.run(
            cmd,
            cwd=IPCR_SOURCE_DIR,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        try:
            tmp_binary.unlink(missing_ok=True)
        except OSError:
            pass
        stderr = getattr(exc, "stderr", "") or ""
        detail = stderr.strip() or str(exc)
        raise IPCRBackendError(
            f"failed to build ipcr from submodule: {detail}"
        ) from exc

    if result.stderr:
        logger.debug("ipcr build stderr: %s", result.stderr.strip())

    os.replace(tmp_binary, IPCR_BINARY)
    if revision is not None:
        IPCR_REVISION_FILE.write_text(revision + "\n", encoding="utf-8")
    else:
        try:
            IPCR_REVISION_FILE.unlink(missing_ok=True)
        except OSError:
            pass
    return IPCR_BINARY


def _resolve_override(candidate: Path) -> Optional[Path]:
    if candidate.is_file():
        return candidate.resolve()
    found = shutil.which(str(candidate))
    return Path(found).resolve() if found else None


def resolve_ipcr_binary(explicit: Optional[Path] = None) -> Path:
    """Resolve ipcr, preferring an explicit override then the pinned submodule."""
    if explicit is not None:
        resolved = _resolve_override(explicit)
        if resolved is None:
            raise IPCRBackendError(f"ipcr executable not found: {explicit}")
        return resolved

    if _source_is_available():
        return build_ipcr_from_submodule()

    # Preserve compatibility with an already-built project-local executable.
    if IPCR_BINARY.is_file():
        return IPCR_BINARY

    system_ipcr = shutil.which("ipcr")
    if system_ipcr:
        return Path(system_ipcr).resolve()

    raise IPCRBackendError(
        "ipcr is unavailable. Initialize the submodule with "
        "'git submodule update --init --recursive', install ipcr on PATH, "
        "or pass --ipcr-bin /path/to/ipcr."
    )
