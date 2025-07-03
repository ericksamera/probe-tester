"""
modules/genome_manager.py

Functions for validating taxa, fetching genome assemblies, and downloading/organizing genome FASTA files.
Supports fetching all child species under a parent taxon (e.g., genus) or a specific species only.
"""

import json
import logging
from typing import Optional, Dict, List, Iterator
from pathlib import Path
import zipfile
from io import TextIOWrapper
from Bio import SeqIO

from modules import io_tools

logger = logging.getLogger(__name__)

def chunked(iterable: List[str], n: int) -> Iterator[List[str]]:
    """Yield successive n-sized chunks from iterable."""
    for i in range(0, len(iterable), n):
        yield iterable[i:i + n]

def validate_taxon(input_taxon: str, parent_mode: bool = False) -> Optional[tuple[str, int, str, dict]]:
    """
    Validates a taxon name using NCBI datasets CLI.
    Returns tuple of (name, taxid, rank, taxonomy dict) for use in both parent/species mode.
    In parent_mode, returns the genus/family/parent (last 'parents' entry).
    """
    logger.info("Validating taxon: %s", input_taxon)
    try:
        result = io_tools.run_command(
            ["datasets", "summary", "taxonomy", "taxon", input_taxon],
            capture_output=True
        )
    except Exception as e:
        logger.error("Error running datasets CLI: %s", e)
        return None

    if not result or not result.stdout:
        logger.warning("No output received for taxon: %s", input_taxon)
        return None

    try:
        result_dict = json.loads(result.stdout.decode())
    except json.JSONDecodeError as e:
        logger.error("Failed to parse JSON: %s", e)
        return None

    reports = result_dict.get("reports")
    if not reports:
        logger.warning("No taxonomy reports found for: %s", input_taxon)
        return None

    taxonomy = reports[0]["taxonomy"]
    current_name = taxonomy["current_scientific_name"]["name"]
    current_taxid = taxonomy["tax_id"]
    rank = taxonomy.get("rank", "").lower()
    parents = taxonomy.get("parents", [])

    # If in parent_mode, use the last parent as the parent taxid (usually the genus)
    if parent_mode and parents:
        parent_taxid = parents[-1]
        classification = taxonomy.get("classification", {})
        # Find the name for this parent taxid in classification
        parent_name = None
        for entry in classification.values():
            if entry["id"] == parent_taxid:
                parent_name = entry["name"]
                break
        # fallback: just call it "parent" if not found
        parent_name = parent_name or f"taxid_{parent_taxid}"
        logger.info("Parent mode: using parent %s (taxid: %s)", parent_name, parent_taxid)
        return parent_name, parent_taxid, "parent", taxonomy

    logger.info("Validated taxon: %s (taxid: %s, rank: %s)", current_name, current_taxid, rank)
    return current_name, current_taxid, rank, taxonomy



def get_species_accessions_by_parent_taxid(
        parent_taxid: int,
        output_dir: Path = Path(".")
    ) -> Optional[Dict[str, List[str]]]:
    """
    Fetch and return mapping of all child species and their accessions under a parent taxon.
    Also writes a CSV with mapping info.
    Returns {species_name: [accession1, accession2, ...]}
    """
    logger.info("Fetching all child species for parent taxid: %d", parent_taxid)
    return get_genomes_mapping(parent_taxid, output_dir)

def get_accessions_for_species(
        species_taxid: int,
        output_dir: Path = Path(".")
    ) -> List[str]:
    """
    Fetch and return accessions for a specific species only.
    Returns [accession1, accession2, ...]
    """
    logger.info("Fetching genome accessions for species taxid: %d", species_taxid)
    mapping = get_genomes_mapping(species_taxid, output_dir)
    if not mapping:
        logger.warning("No genomes found for species taxid: %d", species_taxid)
        return []

    for accessions in mapping.values():
        return accessions
    return []

def get_genomes_mapping(
        taxid: int,
        output_dir: Path = Path(".")
    ) -> Optional[Dict[str, List[str]]]:
    """
    Returns a mapping from species to a list of current_accession IDs for available genomes,
    and writes a CSV file of all assemblies found.
    """
    try:
        result = io_tools.run_command([
                "datasets", "summary", "genome", "taxon", str(taxid),
                "--assembly-source", "GenBank",
                "--assembly-version", "latest",
            ],
            capture_output=True
        )
    except Exception as e:
        logger.error("Error running datasets CLI: %s", e)
        return None

    if not result or not result.stdout:
        logger.warning("No genome assemblies found for taxid: %s", taxid)
        return None

    entries = json.loads(result.stdout.strip())
    sample_mapping = []
    sample_species_dict = {}
    for entry in entries['reports']:
        organism_name = "-".join(entry.get("organism", {}).get("organism_name", "").split()[:2])
        assembly_name: str = entry.get("assembly_info", {}).get("assembly_name")

        if not assembly_name.startswith("ASM"):
            continue

        if organism_name not in sample_species_dict:
            sample_species_dict[organism_name] = []

        sample_species_dict[organism_name].append(entry.get("current_accession"))
        sample_mapping.append({
            "organism_name": organism_name,
            "current_accession": entry.get("current_accession"),
            "assembly_name": assembly_name
        })

    return sample_species_dict

def download_genomes(
        species: str,
        accession_list: List[str],
        work_dir: Path = Path('.'),
        chunk_size: int = 30
    ) -> int:
    """
    Downloads and saves genome FASTA files for the given accession list,
    downloading in chunks of chunk_size if needed.
    Returns the number of genome FASTA files written.
    """
    genomes_dir = work_dir / "genomes" / species
    genomes_dir.mkdir(parents=True, exist_ok=True)

    total = len(accession_list)
    logger.info("Downloading %d genomes for species: %s (chunk size %d)", total, species, chunk_size)
    genomes_counter = 0

    if total <= chunk_size:
        accession_chunks = [accession_list]
    else:
        accession_chunks = list(chunked(accession_list, chunk_size))

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
            bar_task = progress.add_task(f"[cyan]Writing {species} FASTAs", total=total)
            for chunk_idx, acc_chunk in enumerate(accession_chunks):
                logger.info("Downloading chunk %d/%d (%d genomes)", chunk_idx+1, len(accession_chunks), len(acc_chunk))
                try:
                    result = io_tools.run_command(
                        ["datasets", "download", "genome", "accession"] + acc_chunk,
                        capture_output=True
                    )
                except Exception as e:
                    logger.error("Error running datasets CLI on chunk: %s", e)
                    continue

                try:
                    with zipfile.ZipFile("ncbi_dataset.zip") as z:
                        fasta_names = [name for name in z.namelist() if name.endswith(".fna")]
                        for name in fasta_names:
                            with z.open(name) as handle:
                                fasta_io = TextIOWrapper(handle)
                                records = []
                                for i, record in enumerate(SeqIO.parse(fasta_io, "fasta")):
                                    old_id = record.id
                                    record.id = f"{species}.{genomes_counter}.{i}"
                                    record.description = old_id
                                    records.append(record)
                                SeqIO.write(records, genomes_dir / f"{species}.{genomes_counter}.fna", "fasta")
                            genomes_counter += 1
                            progress.update(bar_task, advance=1)
                except Exception as e:
                    logger.error("Error processing FASTA files in chunk: %s", e)
                    continue
    else:
        print(f"[INFO] Writing {species} FASTAs ({total} total):")
        for chunk_idx, acc_chunk in enumerate(accession_chunks):
            print(f"  Downloading chunk {chunk_idx+1}/{len(accession_chunks)} ({len(acc_chunk)} genomes)")
            try:
                result = io_tools.run_command(
                    ["datasets", "download", "genome", "accession"] + acc_chunk,
                    capture_output=True
                )
            except Exception as e:
                print(f"[ERROR] Error running datasets CLI on chunk: {e}")
                continue

            try:
                with zipfile.ZipFile("ncbi_dataset.zip") as z:
                    fasta_names = [name for name in z.namelist() if name.endswith(".fna")]
                    for name in fasta_names:
                        with z.open(name) as handle:
                            accession: str = Path(name).stem.split('_')[0]
                            fasta_io = TextIOWrapper(handle)
                            records = []
                            for i, record in enumerate(SeqIO.parse(fasta_io, "fasta")):
                                old_id = record.id
                                record.id = f"{accession}.{i}"
                                record.description = old_id
                                records.append(record)
                            SeqIO.write(records, genomes_dir / f"{accession}.fna", "fasta")
                        genomes_counter += 1
                        print(f"    [{genomes_counter}/{total}] {name}")
            except Exception as e:
                print(f"[ERROR] Error processing FASTA files in chunk: {e}")
                continue
    
    Path("ncbi_dataset.zip").unlink()
    logger.info("Wrote %d FASTA files for %s", genomes_counter, species)
    return genomes_counter


# ---