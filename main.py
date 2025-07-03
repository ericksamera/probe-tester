# main.py
from pathlib import Path
import logging
import argparse

from modules import genome_downloader, io_tools

def setup_logging(outdir: Path, verbose: bool = False) -> None:
    """Configure logging to both console and a plain log file."""
    log_path = outdir / "pipeline.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    # Console handler
    console_level = logging.DEBUG if verbose else logging.INFO
    c_handler = logging.StreamHandler()
    c_handler.setLevel(console_level)
    c_handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    root.addHandler(c_handler)

    # File handler
    f_handler = logging.FileHandler(log_path, encoding="utf-8")
    f_handler.setLevel(logging.DEBUG)
    f_handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"))
    root.addHandler(f_handler)

    logging.getLogger(__name__).debug("Logging initialised → %s", log_path)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run pipeline with optional command-line input."
    )
    parser.add_argument(
        "-o", "--outdir",
        type=Path,
        default=Path('.'),
        help="Output directory for logs and results"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging to console"
    )
    parser.add_argument(
        "-m", "--max-genomes",
        type=int,
        default=10,
        help="Output"
    )
    parser.add_argument(
        "-d", "--dry-run",
        action="store_true",
        help="Dry run"
    )
    return parser.parse_args()

def main() -> None:
    """
    """

    args = parse_args()

    setup_logging(args.outdir, verbose=args.verbose)

    target_taxon: str = io_tools.ask("Please enter a target taxon.", default="Mycoplasmopsis bovis")

    validated_taxon = genome_downloader.validate_taxon(target_taxon)
    if validated_taxon:
        genomes_mapping = genome_downloader.get_genomes_mapping(validated_taxon[-1])
        if genomes_mapping:
            for species, genomes_list in genomes_mapping.items():
                genome_downloader.download_genomes(species, genomes_list)

    return None

if __name__=="__main__":
    main()
# ---