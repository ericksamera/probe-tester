# modules/fasta_io.py
"""
Minimal FASTA reader/writer with gzip support (no external deps).

API:
  - read_fasta(path, keep_description=True) -> iterator of (id, desc, seq)
  - write_fasta(records, path=None, width=80)
  - write_record(handle, rid, seq, desc=None, width=80)
  - open_guess(path, mode) -> TextIO/BinaryIO for plain/gz/'-'
"""
from __future__ import annotations

from gzip import GzipFile
import gzip
import io
import os
import sys
from typing import (
    Iterable,
    Iterator,
    Literal,
    Optional,
    Tuple,
    Union,
    overload,
    TextIO,
    BinaryIO,
    cast,
)

PathLike = Union[str, os.PathLike]
FileOrPath = Union[PathLike, TextIO]

__all__ = [
    "open_guess",
    "read_fasta",
    "write_fasta",
    "write_record",
]


@overload
def open_guess(
    path: Optional[PathLike],
    mode: Literal["rt", "wt"],
    encoding: str = "utf-8",
) -> TextIO: ...
@overload
def open_guess(
    path: Optional[PathLike],
    mode: Literal["rb", "wb"],
    encoding: str = "utf-8",
) -> BinaryIO: ...
def open_guess(
    path: Optional[PathLike],
    mode: str = "rt",
    encoding: str = "utf-8",
) -> TextIO | BinaryIO:
    """
    Open plain/gz path or '-' for stdio.

    Supports text ('rt'/'wt') and binary ('rb'/'wb') modes. In text mode we
    always return a TextIOWrapper (typed as TextIO).
    """
    is_read = "r" in mode
    is_text = "t" in mode

    # stdio
    if path in (None, "-", ""):
        if is_text:
            buf = sys.stdin.buffer if is_read else sys.stdout.buffer
            return io.TextIOWrapper(buf, encoding=encoding)
        return sys.stdin.buffer if is_read else sys.stdout.buffer

    p = os.fspath(path)

    # gzip
    if p.endswith(".gz"):
        if is_text:
            gz = cast(GzipFile, gzip.open(p, "rb" if is_read else "wb"))
            return io.TextIOWrapper(gz, encoding=encoding)
        return cast(BinaryIO, gzip.open(p, "rb" if is_read else "wb"))

    # plain
    if is_text:
        return cast(TextIO, open(p, mode, encoding=encoding))
    return cast(BinaryIO, open(p, mode))


def read_fasta(
    path: FileOrPath,
    keep_description: bool = True,
) -> Iterator[Tuple[str, str, str]]:
    """
    Stream-parse FASTA from a file path (.gz ok) OR an open text handle.
    Yields (id, desc, seq). If keep_description=False, desc == id.
    """
    # If given an open text handle, use it directly (don't close it here)
    if hasattr(path, "read"):  # TextIO
        fh = cast(TextIO, path)
        _close = False
    else:
        fh = open_guess(cast(PathLike, path), "rt")
        _close = True

    try:
        header: Optional[str] = None
        seq_chunks: list[str] = []
        for line in fh:
            if not line:
                continue
            ch = line[0]
            if ch == ">":
                if header is not None:
                    seq = _join_seq(seq_chunks)
                    rid, desc = _split_header(header, keep_description=keep_description)
                    yield rid, desc, seq
                header = line[1:].strip()
                seq_chunks = []
            elif ch == ";" and header is None:
                continue
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            seq = _join_seq(seq_chunks)
            rid, desc = _split_header(header, keep_description=keep_description)
            yield rid, desc, seq
    finally:
        if _close:
            fh.close()

def write_fasta(
    records: Iterable[Union[Tuple[str, str], Tuple[str, str, str]]],
    path: Optional[PathLike] = None,
    width: int = 80,
) -> None:
    """
    Write records to FASTA. Each record may be (id, seq) or (id, desc, seq).
    If 2-tuple, desc=id. If path is None or '-', writes to stdout.
    """
    with open_guess(path, "wt") as out:
        for rec in records:
            if len(rec) == 2:  # type: ignore[truthy-bool]
                rid, seq = rec  # type: ignore[misc]
                desc = rid
            else:
                rid, desc, seq = rec  # type: ignore[misc]
                if not desc:
                    desc = rid
            out.write(f">{desc}\n")
            _write_wrapped(out, (seq or "").strip().upper(), width)


def write_record(
    handle: TextIO,
    rid: str,
    seq: str,
    desc: Optional[str] = None,
    width: int = 80,
) -> None:
    """Write a single FASTA record to an open text handle."""
    d = desc or rid
    handle.write(f">{d}\n")
    _write_wrapped(handle, (seq or "").strip().upper(), width)


# ---------------------- internal helpers ----------------------


def _split_header(header: str, *, keep_description: bool) -> Tuple[str, str]:
    # ID is the first whitespace-delimited token; description can be the full
    # header if requested.
    h = header.strip()
    if not h:
        return "", ""
    parts = h.split(None, 1)
    rid = parts[0]
    desc = h if keep_description else rid
    return rid, desc


def _write_wrapped(out: TextIO, seq: str, width: int) -> None:
    if width > 0:
        for i in range(0, len(seq), width):
            out.write(seq[i : i + width] + "\n")
    else:
        out.write(seq + "\n")


def _join_seq(chunks: list[str]) -> str:
    # Remove whitespace/newlines and uppercase.
    if not chunks:
        return ""
    s = "".join(chunks)
    # spaces/tabs were stripped when appending; keep for completeness
    return s.replace("\r", "").replace("\n", "").upper()
