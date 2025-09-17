"""
modules/probe_analysis.py
"""
import logging
from typing import List, Tuple, Optional
from pathlib import Path
from tempfile import NamedTemporaryFile

from modules import io_tools
from modules.sequence_utils import reverse_complement, count_mismatches

logger = logging.getLogger(__name__)

def process_genome_job(
    args: tuple
) -> tuple[str, str, list]:
    """Worker function for multiprocessing."""
    (
        forward_primer, reverse_primer, probe, primers_file_path,
        species_name, genome_path, mismatch, dry_run
    ) = args

    # Local imports inside worker to be safe with multiprocessing on some platforms
    from modules.sequence_utils import reverse_complement, count_mismatches
    from modules.probe_analysis import match_probe, run_ipcress

    genome_id = genome_path.stem
    products = run_ipcress(primers_file_path, genome_path, mismatches=mismatch, dry_run=dry_run)
    product_results = []
    if products:
        for product in products:
            f_mismatches = count_mismatches(forward_primer, product[:len(forward_primer)])
            r_mismatches = count_mismatches(
                reverse_complement(reverse_primer), product[-len(reverse_primer):]
            )
            if probe:
                probe_mismatches, probe_pos = match_probe(probe, product)
            else:
                probe_mismatches, probe_pos = None, None
            product_results.append({
                "product": product,
                "forward_mismatches": f_mismatches,
                "reverse_mismatches": r_mismatches,
                "probe_mismatches": probe_mismatches,
                "probe_position": probe_pos,
            })
    return (species_name, genome_id, product_results)

def match_probe(probe: str, sequence: str) -> Tuple[int, int]:
    """
    Find the best match (minimum mismatches) for a probe in a target sequence.
    Returns a tuple of (min_mismatches, best_position).
    """
    min_mismatches = len(probe)
    best_pos = -1
    for i in range(len(sequence) - len(probe) + 1):
        window = sequence[i:i+len(probe)]
        mismatches = count_mismatches(probe, window)
        if mismatches < min_mismatches:
            min_mismatches = mismatches
            best_pos = i
    return min_mismatches, best_pos

def run_ipcress(
        primers_path: Path,
        genome_path: Path,
        mismatches: int = 3,
        dry_run: bool = False
    ) -> Optional[List[str]]:
    """
    Run ipcress with the given primers and genome, returning a list of product sequences.
    """
    result = io_tools.run_command([
        "ipcress", str(primers_path), str(genome_path),
        "--mismatch", str(mismatches),
        "--seed", "0",
        "--pretty", "false",
        "--products", "true"
    ], capture_output=True, dry_run=dry_run)
    if not result:
        logger.warning("No ipcress results for %s vs %s", primers_path, genome_path)
        return None
    raw_lines = [
        line for line in result.stdout.decode().splitlines()
        if "ipcress" not in line
    ]

    sequences = []
    i = 0
    while i < len(raw_lines):
        if raw_lines[i].startswith(">"):
            seq_lines = []
            i += 1
            while i < len(raw_lines) and not (raw_lines[i].startswith(">") or raw_lines[i].startswith("--")):
                seq_lines.append(raw_lines[i].strip())
                i += 1
            sequences.append("".join(seq_lines))
        else:
            i += 1

    logger.debug("Parsed %d ipcress products for %s", len(sequences), genome_path)
    return sequences

def write_ipcress_primers_file(
    forward: str,
    reverse: str,
    primer_min: int = 60,
    primer_max: int = 200
) -> Path:
    """
    Writes a temporary IPCRESS-compatible primers file.
    Format used: NAME FORWARD_PRIMER REVERSE_PRIMER MIN MAX
    (Probe matching, if desired, is handled post hoc in Python.)
    """
    temp = NamedTemporaryFile("w+", delete=False, suffix=".primers")
    temp.write(f"PROBE {forward} {reverse} {primer_min} {primer_max}\n")
    temp.flush()
    return Path(temp.name)

def analyze_genome_products(
    forward_primer: str,
    reverse_primer: str,
    probe: Optional[str],
    genomes_dir: Path,
    primer_min: int = 1,
    primer_max: int = 1500,
    mismatch: int = 3,
    dry_run: bool = False,
    threads: int = 1
) -> dict:
    """
    For each genome in genomes_dir, run ipcress and analyze matches.
    Uses multiprocessing if threads > 1.
    Probe is optional; when absent, probe metrics are omitted (None).
    """
    from modules.probe_analysis import write_ipcress_primers_file

    primers_file = write_ipcress_primers_file(
        forward_primer, reverse_primer, primer_min, primer_max
    )

    all_results = {}

    species_dirs = [d for d in genomes_dir.iterdir() if d.is_dir()]
    jobs = []
    for species_dir in species_dirs:
        species_name = species_dir.name
        genome_paths = sorted(species_dir.glob('*.fna'))
        all_results[species_name] = {}
        for genome_path in genome_paths:
            jobs.append((
                forward_primer, reverse_primer, probe, primers_file,
                species_name, genome_path, mismatch, dry_run
            ))

    if threads > 1:
        try:
            import rich
            rich_available = True
        except ImportError:
            rich_available = False

        if rich_available:
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
                bar_task = progress.add_task(f"[cyan]Analyzing genomes", total=len(jobs))
                from concurrent.futures import ProcessPoolExecutor, as_completed
                with ProcessPoolExecutor(max_workers=threads) as executor:
                    future_to_job = {executor.submit(process_genome_job, job): job for job in jobs}
                    for future in as_completed(future_to_job):
                        species_name, genome_id, product_results = future.result()
                        all_results[species_name][genome_id] = product_results
                        progress.update(bar_task, advance=1)
        else:
            # fallback: simple print loop with multiprocessing
            from concurrent.futures import ProcessPoolExecutor, as_completed
            print(f"Analyzing {len(jobs)} genomes with {threads} processes (rich not installed, no bar)")
            with ProcessPoolExecutor(max_workers=threads) as executor:
                for i, future in enumerate(as_completed([executor.submit(process_genome_job, job) for job in jobs]), 1):
                    species_name, genome_id, product_results = future.result()
                    all_results[species_name][genome_id] = product_results
                    print(f"[{i}/{len(jobs)}] {species_name}:{genome_id}")
    else:
        # Serial fallback
        for job in jobs:
            species_name, genome_id, product_results = process_genome_job(job)
            all_results[species_name][genome_id] = product_results

    try:
        primers_file.unlink()
    except Exception:
        pass

    return all_results
