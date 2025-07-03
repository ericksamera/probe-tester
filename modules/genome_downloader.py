# modules/genome_downloader.py
import json
import logging
from typing import Optional
from modules import io_tools
from pathlib import Path
import pandas as pd
import zipfile

from Bio import SeqIO
from io import TextIOWrapper

logger = logging.getLogger(__name__)

def validate_taxon(input_taxon: str) -> Optional[tuple[str, int]]:
    """
    Validates a taxon name using NCBI datasets CLI.

    If found, prints the current scientific name and taxonomy lineage.
    Returns the validated scientific name, or None if not found.
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
    classification = taxonomy.get("classification", {})
    parents = taxonomy.get("parents", [])

    taxid_to_name = {
        value["id"]: value["name"]
        for level, value in classification.items()
    }

    lineage_parts = [taxid_to_name[tid] for tid in parents if tid in taxid_to_name]
    taxonomy_lineage = "; ".join(lineage_parts)

    print(f"\nScientific name: {current_name}")
    print(f"Taxonomy lineage: {taxonomy_lineage}")

    logger.info("Validated: %s → %s", input_taxon, current_name)

    return current_name, parents[-1]

def get_genomes_mapping(taxid: int, output_dir: Path=Path(".")) -> dict | None:
    """
    """

    try:
        result = io_tools.run_command(
            ["datasets", "summary", "genome", "taxon", str(taxid)],
            capture_output=True
        )
    except Exception as e:
        logger.error("Error running datasets CLI: %s", e)
        return None
    
    if not result or not result.stdout:
        #logger.warning("No output received for taxon: %s", input_taxon)
        return None

    entries: dict = json.loads(result.stdout.strip())
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

    pd.DataFrame(sample_mapping).to_csv(output_dir / ".genome_mapping.csv", index=False)

    return sample_species_dict

def download_genomes(species: str, species_list: list[str], work_dir: Path = Path('.')):
    """
    Downloads and parses genome FASTA files from NCBI and saves standardized .fna files.
    """
    genomes_dir = work_dir / "genomes" / species
    genomes_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading {len(species_list)} genomes for {species} ...")

    try:
        result = io_tools.run_command(
            ["datasets", "download", "genome", "accession"] + species_list,
            capture_output=True
        )
    except Exception as e:
        logger.error("Error running datasets CLI: %s", e)
        return None

    genomes_counter = 0
    try:
        with zipfile.ZipFile("ncbi_dataset.zip") as z:
            for name in z.namelist():
                if name.endswith(".fna"):
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

    except Exception as e:
        logger.error("Error processing FASTA files: %s", e)
        return None

    return genomes_counter
# ---