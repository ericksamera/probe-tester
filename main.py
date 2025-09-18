#!/usr/bin/env python3

__VERSION__ = "1.3.0"
__DESCRIPTION__ = "Functional pipeline for testing amplification."
__AUTHOR__ = "Erick Samera (erick.samera@kpu.ca)"

import sys
import platform

import argparse
import logging
import json
from pathlib import Path

from typing import Optional

from modules.formatting import format_table
from modules.genome_manager import (
    validate_taxon, get_species_accessions_by_parent_taxid,
    get_accessions_for_species, download_genomes
)
from modules.probe_analysis import analyze_genome_products


def check_dependencies(engine: Optional[str] = None, ipcr_bin: Optional[Path] = None, ipcress_bin: Optional[Path] = None,
                       required=["rich"], required_cli=["datasets"]) -> None:
    """
    Verify Python and CLI dependencies. If `engine` is provided, also checks the chosen assay engine executable.
    """
    missing = []
    for mod in required:
        try:
            __import__(mod)
        except ImportError:
            missing.append(mod)

    # Always required for NCBI downloads (datasets)
    from shutil import which
    for cli in required_cli:
        if which(cli) is None:
            missing.append(f"{cli} (CLI)")

    # Engine-specific checks (only when running `assay`)
    if engine:
        if engine == "ipcr":
            exe = str(ipcr_bin) if ipcr_bin else "ipcr"
            if which(exe) is None:
                # try project-local bin/ipcr
                local = Path(__file__).resolve().parent / "bin" / ("ipcr.exe" if sys.platform.startswith("win") else "ipcr")
                if not local.exists():
                    missing.append("ipcr (CLI)")
        elif engine == "ipcress":
            exe = str(ipcress_bin) if ipcress_bin else "ipcress"
            if which(exe) is None:
                missing.append("ipcress (CLI)")

    if missing:
        print("[ERROR] Missing dependencies:")
        for item in missing:
            print(f"  - {item}")
        print("Please install required Python modules and CLI tools before running.")
        sys.exit(1)


def list_command(args):
    parent_mode = args.mode == "parent"
    res = validate_taxon(args.taxon, parent_mode=parent_mode)
    if not res:
        logging.error("Could not validate taxon: %s", args.taxon)
        return
    sci_name, taxid, rank, taxonomy = res

    if args.mode == "parent":
        mapping = get_species_accessions_by_parent_taxid(taxid)
        if not mapping:
            logging.error("No genomes found under parent taxon.")
            return
        print(f"\nSpecies under parent {sci_name} (taxid {taxid}):")
        for species, accessions in mapping.items():
            print(f"  {species}: {len(accessions)} genomes")
    else:
        accessions = get_accessions_for_species(taxid)
        print(f"\nSpecies: {sci_name} (taxid {taxid})")
        print(f"  {len(accessions)} genomes")
        for acc in accessions:
            print(f"    {acc}")


def download_command(args):
    parent_mode = args.mode == "parent"
    res = validate_taxon(args.taxon, parent_mode=parent_mode)
    if not res:
        logging.error("Could not validate taxon: %s", args.taxon)
        return
    sci_name, taxid, rank, taxonomy = res

    if args.mode == "parent":
        mapping = get_species_accessions_by_parent_taxid(taxid, output_dir=args.outdir)
        if not mapping:
            logging.error("No genomes found under parent taxon.")
            return
        species_list = list(mapping.items())
        logging.info("Would download the following species and genome counts:")
        for species, accessions in species_list:
            count = len(accessions)
            limited = f" (limited to {args.max_genomes})" if args.max_genomes and count > args.max_genomes else ""
            print(f"  {species}: {min(count, args.max_genomes or count)} genomes{limited}")

        # Prompt for confirmation if a large download
        species_total = len(species_list)
        genome_total = sum(len(accessions) if not args.max_genomes else min(len(accessions), args.max_genomes)
                           for _, accessions in species_list)
        BIG_DOWNLOAD_THRESHOLD = 500
        if not getattr(args, "force", False) and genome_total >= BIG_DOWNLOAD_THRESHOLD:
            print(f"[WARNING] You are about to download {genome_total} genomes ({species_total} species).")
            ans = input("Proceed? [y/N]: ").strip().lower()
            if not (ans == "y" or ans == "yes"):
                print("Aborted.")
                return

        if not args.dry_run:
            try:
                import rich
                rich_available = True
            except ImportError:
                rich_available = False

            if rich_available:
                from rich.progress import Progress, BarColumn, TimeElapsedColumn, TimeRemainingColumn, TextColumn
                from rich.console import Console
                total_fastas = sum(min(len(acc), args.max_genomes) if args.max_genomes else len(acc)
                                  for _, acc in species_list)
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
                    overall_task = progress.add_task("[magenta]Total FASTAs", total=total_fastas)
                    for species, accessions in species_list:
                        if args.max_genomes:
                            accessions = accessions[:args.max_genomes]
                        download_genomes(
                            species,
                            accessions,
                            work_dir=args.outdir,
                            progress=progress,
                            overall_task=overall_task
                        )
                        logging.debug("Downloaded %d genomes for %s", len(accessions), species)
            else:
                print(f"[INFO] Downloading {len(species_list)} species:")
                for i, (species, accessions) in enumerate(species_list, 1):
                    if args.max_genomes:
                        accessions = accessions[:args.max_genomes]
                    print(f"  [{i}/{len(species_list)}] {species} ({len(accessions)} genomes)")
                    download_genomes(species, accessions, work_dir=args.outdir)
                    logging.info("Downloaded %d genomes for %s", len(accessions), species)
        print(f"\n[SUCCESS] Run complete!\nNext: \n python main.py assay --forward <FWD> --reverse <REV> --probe <PROBE>")
    else:  # mode == "species"
        accessions = get_accessions_for_species(taxid, output_dir=args.outdir)
        if not accessions:
            logging.error("No genomes found for the given species.")
            return
        count = len(accessions)
        limited = f" (limited to {args.max_genomes})" if args.max_genomes and count > args.max_genomes else ""
        print(f"Would download: {min(count, args.max_genomes or count)} genomes for {sci_name}{limited}")
        if not args.dry_run:
            if args.max_genomes:
                accessions = accessions[:args.max_genomes]
            download_genomes(sci_name.replace(" ", "-"), accessions, work_dir=args.outdir)
            logging.info("Downloaded %d genomes for %s", len(accessions), sci_name)
        print(f"\n[SUCCESS] Run complete!\nNext: \n python main.py assay --forward <FWD> --reverse <REV> --probe <PROBE>")


def assay_command(args):
    import datetime
    import re

    if args.run_name:
        run_name = args.run_name
    else:
        now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        run_name = f"run_{now}"

    slug = re.sub(r'[^A-Za-z0-9._-]', '_', run_name.strip())
    if args.output == Path("results.json"):
        output_path = Path(f"results-{slug}.json")
    else:
        output_path = args.output

    threads = getattr(args, "threads", 1)
    results = analyze_genome_products(
        forward_primer=args.forward,
        reverse_primer=args.reverse,
        probe=args.probe,  # Optional
        genomes_dir=args.genomes_dir,
        primer_min=args.primer_min,
        primer_max=args.primer_max,
        mismatch=args.mismatch,
        threads=threads,
        engine=args.engine,
        ipcr_bin=args.ipcr_bin,
        ipcress_bin=args.ipcress_bin,
    )

    metadata = {
        "run_name": run_name,
        "command": " ".join(sys.argv),
        "date": datetime.datetime.now().isoformat(),
        "forward_primer": args.forward,
        "reverse_primer": args.reverse,
        "probe": args.probe,  # May be None
        "primer_min": args.primer_min,
        "primer_max": args.primer_max,
        "mismatch": args.mismatch,
        "genomes_dir": str(args.genomes_dir),
        "probe_tester_version": __VERSION__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "engine": args.engine,
    }

    output_data = {
        "_metadata": metadata,
        run_name: results
    }

    with open(output_path, "w") as f:
        json.dump(output_data, f, indent=2)

    logging.info("Results written to %s under run name '%s'", output_path, run_name)
    print(f"\n[SUCCESS] Run complete!\nResults written to {output_path}\nNext: python main.py summarize --input {output_path}")


def summarize_command(args):
    summarize_results(args.input, getattr(args, "format", "text"), args.target)


def setup_logging(outdir: Path, verbose: bool = False) -> None:
    """Configure logging to both console and a plain log file."""
    log_path = outdir / "pipeline.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    root = logging.getLogger()
    root.setLevel(logging.DEBUG if verbose else logging.INFO)

    # Console handler
    c_handler = logging.StreamHandler()
    c_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    c_handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    root.addHandler(c_handler)

    # File handler
    f_handler = logging.FileHandler(log_path, encoding="utf-8")
    f_handler.setLevel(logging.DEBUG)
    f_handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"))
    root.addHandler(f_handler)

    logging.getLogger(__name__).debug("Logging initialized at %s", log_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"probe-tester v{__VERSION__} | {__DESCRIPTION__}",
        epilog=f"{__AUTHOR__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- LIST SUBCOMMAND ---
    listcmd = subparsers.add_parser("list", help="List available species/genomes for a taxon")
    listcmd.add_argument("--taxon", required=True, help="NCBI taxon name or taxid (parent or species)")
    listcmd.add_argument("--mode", choices=["parent", "species"], default="parent",
                         help="List all genomes under a parent taxon or just a single species")

    # --- DOWNLOAD SUBCOMMAND ---
    dl = subparsers.add_parser("download", help="Download genomes")
    dl.add_argument("--taxon", required=True, help="NCBI taxon name or taxid (parent or species)")
    dl.add_argument("--mode", choices=["parent", "species"], default="parent",
                    help="Download all genomes under a parent taxon or just a single species")
    dl.add_argument("--outdir", type=Path, default=Path("."), help="Output directory")
    dl.add_argument("--max-genomes", type=int, default=None, help="Limit number of genomes per species (optional)")
    dl.add_argument("--verbose", action="store_true")
    dl.add_argument("--dry-run", action="store_true", help="Preview what would be downloaded, but do not download")

    # --- TEST SUBCOMMAND ---
    test = subparsers.add_parser("assay", help="Test probes on downloaded genomes")
    test.add_argument("--forward", required=True, help="Forward primer sequence")
    test.add_argument("--reverse", required=True, help="Reverse primer sequence")
    test.add_argument("--probe", required=False, help="Probe sequence (optional)")
    test.add_argument("--genomes-dir", type=Path, default=Path("./genomes"), help="Genomes directory (output of download)")
    test.add_argument("--primer-min", type=int, default=60, help="Min product length")
    test.add_argument("--primer-max", type=int, default=200, help="Max product length")
    test.add_argument("--mismatch", type=int, default=3, help="Max mismatches per primer")
    test.add_argument("--output", type=Path, default=Path("results.json"), help="Output JSON results file")
    test.add_argument("--run-name", type=str, default=None, help="Name for this set of primers/probe results (optional)")
    test.add_argument("--verbose", action="store_true")
    test.add_argument("--threads", type=int, default=1, help="Number of processes to use for parallel genome analysis")

    # Engine selection
    test.add_argument("--engine", choices=["ipcr", "ipcress"], default="ipcr",
                      help="Which assay engine to use (default: ipcr)")
    test.add_argument("--ipcr-bin", type=Path, default=None,
                      help="Path to ipcr executable (if not in ./bin or PATH)")
    test.add_argument("--ipcress-bin", type=Path, default=None,
                      help="Path to ipcress executable (optional)")

    # --- SUMMARIZE SUBCOMMAND ---
    summ = subparsers.add_parser("summarize", help="Summarize probe test results")
    summ.add_argument("--input", type=Path, default=Path("results.json"), help="Input JSON results file")
    summ.add_argument("--format", choices=["text", "csv", "markdown"], default="text", help="Output format")
    summ.add_argument("--target", nargs="+", help="Organism(s) to treat as targets (exact or glob match, e.g. 'Mycoplasmopsis-bovis' or 'Mycoplasmopsis-*')")

    return parser.parse_args()


def summarize_results(results_path: Path, output_format: str = "text", target: Optional[list] = None) -> None:
    """Print detailed summary per organism from a results JSON, in text/csv/markdown, with totals and proper case-insensitive panel separation."""
    import csv
    from io import StringIO
    import fnmatch

    with results_path.open() as input_file:
        data = json.load(input_file)

    run_items = [(k, v) for k, v in data.items() if not k.startswith("_")]

    def build_rows_and_totals(run_data, selected_organisms):
        rows = []
        grand_fwd = grand_rev = grand_probe = grand_amplicons = grand_genomes = grand_genomes_with_hits = 0
        grand_probe_n = 0
        for organism in selected_organisms:
            genomes = run_data[organism]
            total_fwd = total_rev = total_probe = total_amplicons = 0
            total_probe_n = 0
            n_genomes = n_genomes_with_hits = 0
            for genome_results in genomes.values():
                n_genomes += 1
                if genome_results:
                    n_genomes_with_hits += 1
                    total_amplicons += len(genome_results)
                    for entry in genome_results:
                        total_fwd += entry.get("forward_mismatches", 0)
                        total_rev += entry.get("reverse_mismatches", 0)
                        pm = entry.get("probe_mismatches", None)
                        if pm is not None:
                            total_probe += pm
                            total_probe_n += 1
            avg_fwd = total_fwd / total_amplicons if total_amplicons else "-"
            avg_rev = total_rev / total_amplicons if total_amplicons else "-"
            avg_probe = (total_probe / total_probe_n) if total_probe_n else "-"
            avg_amplicons_per_genome = (
                total_amplicons / n_genomes_with_hits
                if n_genomes_with_hits else "-"
            )
            fraction_with_hits = (n_genomes_with_hits / n_genomes * 100) if n_genomes else 0

            rows.append([
                organism,
                f"{avg_fwd:.2f}" if isinstance(avg_fwd, float) else "-",
                f"{avg_rev:.2f}" if isinstance(avg_rev, float) else "-",
                f"{avg_probe:.2f}" if isinstance(avg_probe, float) else "-",
                str(total_amplicons),
                str(n_genomes),
                f"{avg_amplicons_per_genome:.2f}" if isinstance(avg_amplicons_per_genome, float) else "-",
                str(n_genomes_with_hits),
                f"{fraction_with_hits:.1f}%"
            ])

            grand_fwd += total_fwd
            grand_rev += total_rev
            grand_probe += total_probe
            grand_amplicons += total_amplicons
            grand_genomes += n_genomes
            grand_genomes_with_hits += n_genomes_with_hits
            grand_probe_n += total_probe_n

        grand_avg_fwd = grand_fwd / grand_amplicons if grand_amplicons else "-"
        grand_avg_rev = grand_rev / grand_amplicons if grand_amplicons else "-"
        grand_avg_probe = (grand_probe / grand_probe_n) if grand_probe_n else "-"
        grand_avg_amp_per_genome = (
            grand_amplicons / grand_genomes_with_hits
            if grand_genomes_with_hits else "-"
        )
        grand_frac_with_hits = (grand_genomes_with_hits / grand_genomes * 100) if grand_genomes else 0

        total_row = [
            "[TOTAL]",
            f"{grand_avg_fwd:.2f}" if isinstance(grand_avg_fwd, float) else "-",
            f"{grand_avg_rev:.2f}" if isinstance(grand_avg_rev, float) else "-",
            f"{grand_avg_probe:.2f}" if isinstance(grand_avg_probe, float) else "-",
            str(grand_amplicons),
            str(grand_genomes),
            f"{grand_avg_amp_per_genome:.2f}" if isinstance(grand_avg_amp_per_genome, float) else "-",
            str(grand_genomes_with_hits),
            f"{grand_frac_with_hits:.1f}%"
        ]
        rows.append(total_row)
        return rows

    headers = [
        "Organism",
        "Fwd MM",
        "Rev MM",
        "Probe MM",
        "Amplicons",
        "Genomes Tested",
        "Avg. Amp./Genome",
        "Genomes w/ Amp.",
        "% Genomes w/ Amp."
    ]

    for run_name, run_data in run_items:
        all_organisms = list(run_data.keys())

        # Map lowercase organism -> canonical
        org_lc_to_canon = {org.lower(): org for org in all_organisms}
        patterns = [p.lower() for p in target] if target else []
        # Find targets case-insensitively, using fnmatch wildcards
        target_set = set()
        for org_lc, org in org_lc_to_canon.items():
            if patterns and any(fnmatch.fnmatchcase(org_lc, pat) for pat in patterns):
                target_set.add(org)
        target_organisms = sorted(target_set)
        nontarget_organisms = sorted([org for org in all_organisms if org not in target_set])

        # Panel 1: targets
        if target_organisms:
            print(f"\n=== Run: {run_name} (Targets) ===")
            target_rows = build_rows_and_totals(run_data, target_organisms)
            if output_format == "csv":
                output = StringIO()
                writer = csv.writer(output)
                writer.writerow(["Run", *headers])
                for row in target_rows:
                    writer.writerow([run_name, *row])
                print(output.getvalue())
            elif output_format == "markdown":
                print(f"#### TARGET ORGANISMS")
                md = "| " + " | ".join(headers) + " |\n"
                md += "|---" * len(headers) + "|\n"
                for row in target_rows:
                    md += "| " + " | ".join(row) + " |\n"
                print(md)
            else:
                print(format_table(headers=headers, rows=target_rows))
        
        # Panel 2: nontargets
        if nontarget_organisms:
            print(f"\n=== Run: {run_name} {'(Non-targets)' if target else ''} ===")
            nontarget_rows = build_rows_and_totals(run_data, nontarget_organisms)
            if output_format == "csv":
                output = StringIO()
                writer = csv.writer(output)
                writer.writerow(["Run", *headers])
                for row in nontarget_rows:
                    writer.writerow([run_name, *row])
                print(output.getvalue())
            elif output_format == "markdown":
                print(f"#### NON-TARGET ORGANISMS")
                md = "| " + " | ".join(headers) + " |\n"
                md += "|---" * len(headers) + "|\n"
                for row in nontarget_rows:
                    md += "| " + " | ".join(row) + " |\n"
                print(md)
            else:
                print(format_table(headers=headers, rows=nontarget_rows))


def main():
    args = parse_args()

    # Engine-aware dependency check (only relevant for `assay`)
    if args.command == "assay":
        check_dependencies(engine=args.engine, ipcr_bin=args.ipcr_bin, ipcress_bin=args.ipcress_bin)
    else:
        check_dependencies()

    setup_logging(args.outdir if hasattr(args, "outdir") else Path("."), verbose=getattr(args, 'verbose', False))

    # Dispatch to modular command functions
    if args.command == "list":
        list_command(args)
    elif args.command == "download":
        download_command(args)
    elif args.command == "assay":
        assay_command(args)
    elif args.command == "summarize":
        summarize_command(args)
    else:
        print("[ERROR] Unknown command.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[INFO] Cancelled by user. Partial results may be saved.")
        sys.exit(1)
