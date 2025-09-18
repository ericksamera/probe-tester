"""
modules/genome_manager.py

Functions for validating taxa, fetching genome assemblies, and downloading/organizing genome FASTA files.
Supports fetching all child species under a parent taxon (e.g., genus) or a specific species only.
"""

import json
import logging
from typing import Any, List, Dict, Iterator, Optional
from pathlib import Path
import zipfile
import io

from modules import io_tools
from modules.fasta_io import read_fasta, write_fasta

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
    Returns a mapping from species to a list of current_accession IDs for available genomes.
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
    sample_species_dict: Dict[str, List[str]] = {}
    for entry in entries['reports']:
        organism_name = "-".join(entry.get("organism", {}).get("organism_name", "").split()[:2])
        assembly_name: str = entry.get("assembly_info", {}).get("assembly_name")

        if not assembly_name or not assembly_name.startswith("ASM"):
            continue

        sample_species_dict.setdefault(organism_name, []).append(entry.get("current_accession"))

    return sample_species_dict

def download_genomes(
    species: str,
    accession_list: List[str],
    work_dir: Path = Path("."),
    chunk_size: int = 30,
    progress: Optional[Any] = None,
    overall_task: Optional[Any] = None,
) -> int:
    """
    Download genome FASTAs for `accession_list` into {work_dir}/genomes/{species}.
    Returns the number of FASTA files written.

    progress/overall_task: optional Rich objects (typed as Any to avoid dependency).
    """
    genomes_dir: Path = work_dir / "genomes" / species
    genomes_dir.mkdir(parents=True, exist_ok=True)
    zip_path: Path = work_dir / "ncbi_dataset.zip"

    total: int = len(accession_list)
    genomes_counter: int = 0

    accession_chunks: List[List[str]] = (
        [accession_list] if total <= chunk_size else list(chunked(accession_list, chunk_size))
    )

    def _extract_accession(member_path_str: str) -> str:
        p = Path(member_path_str)
        candidates = [p.parent.name, p.stem, p.stem.split("_genomic")[0]]
        for c in candidates:
            if c.startswith(("GCA_", "GCF_")):
                return c
        return candidates[-1]

    def _process_zip() -> None:
        nonlocal genomes_counter
        with zipfile.ZipFile(zip_path) as z:
            fasta_names: List[str] = [name for name in z.namelist() if name.endswith(".fna")]
            for name in fasta_names:
                with z.open(name) as handle, io.TextIOWrapper(handle) as fasta_text:
                    accession: str = _extract_accession(name)
                    records: List[tuple[str, str, str]] = []
                    for i, (rid, desc, seq) in enumerate(read_fasta(fasta_text)):
                        new_id = f"{accession}.{i}"
                        records.append((new_id, rid, seq))  # (id, desc, seq)
                    write_fasta(records, genomes_dir / f"{accession}.fna")
                genomes_counter += 1

    if progress is not None and overall_task is not None:
        species_task: Any = progress.add_task(f"[cyan]{species} FASTAs", total=total)
        for acc_chunk in accession_chunks:
            try:
                io_tools.run_command(
                    ["datasets", "download", "genome", "accession", *acc_chunk],
                    capture_output=True,
                    cwd=work_dir,
                )
            except Exception as e:
                logger.error("Error running datasets CLI on chunk: %s", e)
                continue

            try:
                _process_zip()
                # advance both bars by count just processed
                progress.update(overall_task, advance=len(acc_chunk))
                progress.update(species_task, advance=len(acc_chunk))
            except Exception as e:
                logger.error("Error processing FASTA files in chunk: %s", e)
            finally:
                try:
                    if zip_path.exists():
                        zip_path.unlink()
                except Exception:
                    pass

        progress.update(species_task, completed=total)
        progress.remove_task(species_task)
    else:
        print(f"[INFO] Writing {species} FASTAs ({total} total):")
        for idx, acc_chunk in enumerate(accession_chunks, 1):
            print(f"  Downloading chunk {idx}/{len(accession_chunks)} ({len(acc_chunk)} genomes)")
            try:
                io_tools.run_command(
                    ["datasets", "download", "genome", "accession", *acc_chunk],
                    capture_output=True,
                    cwd=work_dir,
                )
            except Exception as e:
                print(f"[ERROR] Error running datasets CLI on chunk: {e}")
                continue

            try:
                before = genomes_counter
                _process_zip()
                processed = genomes_counter - before
                for _ in range(processed):
                    print(f"    [{genomes_counter}/{total}] wrote FASTA")
            except Exception as e:
                print(f"[ERROR] Error processing FASTA files in chunk: {e}")
            finally:
                try:
                    if zip_path.exists():
                        zip_path.unlink()
                except Exception:
                    pass

    return genomes_counter