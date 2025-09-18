"""
modules/probe_analysis.py
Drop-in: supports --engine ipcress | ipcr while keeping probe handling in Python.
"""
from __future__ import annotations

import logging
from typing import List, Tuple, Optional, Dict
from pathlib import Path
from tempfile import NamedTemporaryFile
from concurrent.futures import ProcessPoolExecutor, as_completed

from modules import io_tools
from modules.sequence_utils import reverse_complement, count_mismatches

logger = logging.getLogger(__name__)

# Recognize common FASTA file patterns inside each species directory
_FASTA_GLOBS = ("*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz", "*.fasta.gz")


# ----------------------------- FASTA parsing -----------------------------

def _parse_fasta_from_text(text: str, drop_prefix_tokens: Optional[List[str]] = None) -> List[str]:
    """
    Very small FASTA parser for CLI outputs we produce (ipcress/ipcr).
    Returns a list of plain product sequences (uppercased, no gaps).
    """
    seqs: List[str] = []
    lines = text.splitlines()
    # Optionally drop any non‑FASTA chatter lines that contain specific tokens
    if drop_prefix_tokens:
        keep = []
        for ln in lines:
            lower = ln.lower()
            if any(tok in lower for tok in drop_prefix_tokens):
                continue
            keep.append(ln)
        lines = keep

    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln.startswith(">"):
            i += 1
            buf: List[str] = []
            while i < len(lines):
                nxt = lines[i].strip()
                if nxt.startswith(">"):
                    break
                if nxt.startswith("--"):  # some tools end blocks with separator lines
                    break
                if nxt:
                    buf.append(nxt)
                i += 1
            if buf:
                seqs.append("".join(buf).upper())
            continue
        i += 1
    return seqs


# ----------------------------- Probe matching -----------------------------

def match_probe(probe: str, sequence: str) -> Tuple[int, int]:
    """
    Find the best match (minimum mismatches) for a probe in a target sequence.
    Returns (min_mismatches, best_position). If no window returned, (-1,-1).
    """
    if not probe or not sequence or len(sequence) < len(probe):
        return -1, -1
    min_mismatches = len(probe)
    best_pos = -1
    for i in range(len(sequence) - len(probe) + 1):
        window = sequence[i:i + len(probe)]
        mismatches = count_mismatches(probe, window)
        if mismatches < min_mismatches:
            min_mismatches = mismatches
            best_pos = i
            if min_mismatches == 0:
                break
    return min_mismatches, best_pos


# ----------------------------- IPCRESS runner -----------------------------

def write_ipcress_primers_file(
    forward: str,
    reverse: str,
    primer_min: int = 60,
    primer_max: int = 200,
) -> Path:
    """
    Write an IPCRESS-compatible primers file: ID FORWARD REVERSE MIN MAX
    Probe handling stays in Python (post hoc).
    """
    tmp = NamedTemporaryFile("w+", delete=False, suffix=".primers")
    tmp.write(f"PROBE {forward} {reverse} {int(primer_min)} {int(primer_max)}\n")
    tmp.flush()
    return Path(tmp.name)


def run_ipcress(
    primers_path: Path,
    genome_path: Path,
    mismatches: int = 3,
    dry_run: bool = False,
    ipcress_exe: Optional[Path] = None,
) -> Optional[List[str]]:
    """
    Run ipcress and return a list of product sequences (FASTA payload only).
    We request product sequences and suppress pretty/tabular noise.
    """
    exe = str(ipcress_exe) if ipcress_exe else "ipcress"
    result = io_tools.run_command([
        exe, str(primers_path), str(genome_path),
        "--mismatch", str(mismatches),
        "--seed", "0",
        "--pretty", "false",
        "--products", "true",
    ], capture_output=True, dry_run=dry_run)
    if not result:
        logger.warning("No ipcress result for %s vs %s", primers_path, genome_path)
        return None

    # Filter out "ipcress:" lines emitted before FASTA, then parse FASTA
    text = result.stdout.decode(errors="replace")
    seqs = _parse_fasta_from_text(text, drop_prefix_tokens=["ipcress:"])
    logger.debug("ipcress parsed %d products for %s", len(seqs), genome_path)
    return seqs


# ------------------------------- IPCR runner -------------------------------

def run_ipcr(
    forward: str,
    reverse: str,
    genome_path: Path,
    *,
    min_len: int,
    max_len: int,
    mismatches: int = 0,
    threads: int = 1,
    dry_run: bool = False,
    ipcr_exe: Optional[Path] = None,
) -> Optional[List[str]]:
    """
    Run Go 'ipcr' in FASTA mode and return product sequences.

    Exit codes:
      0 -> parse FASTA, return sequences
      1 -> no amplicons found => return []
      2 -> usage/config error => ValueError
      3+ -> runtime/IO error   => RuntimeError
    """
    import subprocess

    exe = str(ipcr_exe) if ipcr_exe else "ipcr"
    cmd = [
        exe,
        "--forward", forward,
        "--reverse", reverse,
        "--sequences", str(genome_path),
        "--min-length", str(int(min_len)),
        "--max-length", str(int(max_len)),
        "--mismatches", str(int(mismatches)),
        "--threads", str(int(threads if threads > 0 else 1)),
        "--output", "fasta",
        "--sort",
    ]

    try:
        result = io_tools.run_command(cmd, capture_output=True, dry_run=dry_run)
        rc = getattr(result, "returncode", 0)
        stdout_b = getattr(result, "stdout", b"") or b""
        stderr_b = getattr(result, "stderr", b"") or b""
    except subprocess.CalledProcessError as e:
        rc = e.returncode
        stdout_b = e.stdout or b""
        stderr_b = e.stderr or b""

    stdout = stdout_b.decode(errors="replace")
    stderr = stderr_b.decode(errors="replace")

    if rc == 0:
        seqs = _parse_fasta_from_text(stdout)
        logger.debug("ipcr parsed %d products for %s", len(seqs), genome_path)
        return seqs
    if rc == 1:
        logger.info("ipcr: no amplicons for %s", genome_path)
        return []
    if rc == 2:
        raise ValueError(
            f"ipcr usage/config error on {genome_path}: {stderr.strip() or 'no stderr'}"
        )
    raise RuntimeError(
        f"ipcr failed on {genome_path} (exit {rc}): {stderr.strip() or 'no stderr'}"
    )

# --------------------------- Per‑genome worker ----------------------------

def process_genome_job(args: tuple) -> tuple[str, str, list]:
    """Worker for multiprocessing."""
    (
        engine, forward_primer, reverse_primer, probe,
        primers_file_path,      # may be None for ipcr
        species_name, genome_path,
        mismatch, primer_min, primer_max,
        dry_run, ipcr_bin, ipcress_bin,
    ) = args

    genome_id = genome_path.stem

    # Run engine → get amplicon sequences
    if engine == "ipcress":
        products = run_ipcress(
            primers_file_path, genome_path,
            mismatches=mismatch, dry_run=dry_run,
            ipcress_exe=ipcress_bin,
        ) or []
    elif engine == "ipcr":
        products = run_ipcr(
            forward_primer, reverse_primer, genome_path,
            min_len=primer_min, max_len=primer_max,
            mismatches=mismatch, threads=1, dry_run=dry_run,
            ipcr_exe=ipcr_bin,
        ) or []
    else:
        raise ValueError(f"unknown engine {engine!r}")

    # Score primers (and optional probe) in Python
    product_results = []
    for prod in products:
        f_mismatches = count_mismatches(forward_primer, prod[:len(forward_primer)])
        r_mismatches = count_mismatches(
            reverse_complement(reverse_primer), prod[-len(reverse_primer):]
        )
        if probe:
            probe_mismatches, probe_pos = match_probe(probe, prod)
        else:
            probe_mismatches, probe_pos = None, None

        product_results.append({
            "product": prod,
            "forward_mismatches": f_mismatches,
            "reverse_mismatches": r_mismatches,
            "probe_mismatches": probe_mismatches,
            "probe_position": probe_pos,
        })

    return (species_name, genome_id, product_results)


# --------------------------- Public entry point ---------------------------

def analyze_genome_products(
    *,
    forward_primer: str,
    reverse_primer: str,
    probe: Optional[str],
    genomes_dir: Path,
    primer_min: int = 60,
    primer_max: int = 200,
    mismatch: int = 3,
    dry_run: bool = False,
    threads: int = 1,
    engine: str = "ipcr",
    ipcr_bin: Optional[Path] = None,
    ipcress_bin: Optional[Path] = None,
) -> Dict[str, Dict[str, list]]:
    """
    For each genome in genomes_dir, run the chosen engine and compute:
      { species: { genome_id: [ {product, forward_mismatches, reverse_mismatches, probe_*}... ] } }

    Notes:
    - Probe handling is done here (not in the engine), matching your original design. :contentReference[oaicite:2]{index=2}
    - For ipcr, we emit FASTA products and parse sequences (Go CLI supports --output fasta). :contentReference[oaicite:3]{index=3}
    """
    engine = (engine or "ipcr").lower()
    if engine not in {"ipcress", "ipcr"}:
        raise ValueError("engine must be 'ipcress' or 'ipcr'")

    # Prepare ipcress primer file if needed
    primers_file: Optional[Path] = None
    if engine == "ipcress":
        primers_file = write_ipcress_primers_file(forward_primer, reverse_primer, primer_min, primer_max)

    # Build job list
    genomes_dir = Path(genomes_dir)
    species_dirs = [d for d in genomes_dir.iterdir() if d.is_dir()]
    jobs = []
    results: Dict[str, Dict[str, list]] = {}

    for species_dir in species_dirs:
        species_name = species_dir.name
        results[species_name] = {}
        genome_paths: List[Path] = []
        for pat in _FASTA_GLOBS:
            genome_paths.extend(species_dir.glob(pat))
        genome_paths = sorted(set(genome_paths))
        for gp in genome_paths:
            jobs.append((
                engine, forward_primer, reverse_primer, probe,
                primers_file,
                species_name, gp,
                mismatch, primer_min, primer_max,
                dry_run, ipcr_bin, ipcress_bin,
            ))

    # Execute
    if threads and threads > 1:
        try:
            import rich  # optional progress UI
            _has_rich = True
        except ImportError:
            _has_rich = False

        if _has_rich:
            from rich.progress import Progress, BarColumn, TimeElapsedColumn, TimeRemainingColumn, TextColumn
            from rich.console import Console
            console = Console()
            with Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("{task.completed}/{task.total}"),
                TimeElapsedColumn(),
                TimeRemainingColumn(),
                refresh_per_second=5,
                console=console,
                transient=False,
            ) as progress:
                bar = progress.add_task("[cyan]Analyzing genomes", total=len(jobs))
                with ProcessPoolExecutor(max_workers=threads) as ex:
                    fut_to_job = {ex.submit(process_genome_job, j): j for j in jobs}
                    for fut in as_completed(fut_to_job):
                        species_name, genome_id, product_results = fut.result()
                        results[species_name][genome_id] = product_results
                        progress.update(bar, advance=1)
        else:
            print(f"Analyzing {len(jobs)} genomes with {threads} processes")
            with ProcessPoolExecutor(max_workers=threads) as ex:
                for i, fut in enumerate(as_completed([ex.submit(process_genome_job, j) for j in jobs]), 1):
                    species_name, genome_id, product_results = fut.result()
                    results[species_name][genome_id] = product_results
                    print(f"[{i}/{len(jobs)}] {species_name}:{genome_id}")
    else:
        # Serial
        for j in jobs:
            species_name, genome_id, product_results = process_genome_job(j)
            results[species_name][genome_id] = product_results

    # Cleanup
    try:
        if primers_file:
            primers_file.unlink(missing_ok=True)
    except Exception:
        pass

    return results
