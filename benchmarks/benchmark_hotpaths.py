#!/usr/bin/env python3
"""Compare these hot paths with an ORIGINAL, TRUSTED codebase JSON export.

Run from the repository root:
    python benchmarks/benchmark_hotpaths.py /path/to/original-codebase.json \
        --output benchmark-results.json

This imports Python code from the supplied export; do not use an untrusted file.
The benchmark uses synthetic sequences, not an ipcr/ipcress executable.
"""

from __future__ import annotations

import argparse
import importlib
import json
import platform
import random
import statistics
import sys
import tempfile
import time
from pathlib import Path


def load_modules(root: Path):
    saved = {
        key: val
        for key, val in sys.modules.items()
        if key == "modules" or key.startswith("modules.")
    }
    for key in saved:
        del sys.modules[key]
    sys.path.insert(0, str(root))
    try:
        su = importlib.import_module("modules.sequence_utils")
        pa = importlib.import_module("modules.probe_analysis")
        return su, pa
    finally:
        sys.path.pop(0)
        for key in list(sys.modules):
            if key == "modules" or key.startswith("modules."):
                del sys.modules[key]
        sys.modules.update(saved)


def time_pair(label, before, after, count, *, repeats=7, loops=1):
    # Warm the relevant caches and require matching results before timing.
    if before() != after():
        raise AssertionError("Result differs: " + label)
    timings = [[], []]
    for repeat in range(repeats):
        for index in ((0, 1) if repeat % 2 == 0 else (1, 0)):
            func = (before, after)[index]
            start = time.perf_counter()
            for _ in range(loops):
                func()
            timings[index].append((time.perf_counter() - start) / loops)
    old, new = map(statistics.median, timings)
    result = {
        "operation": label,
        "items_per_batch": count,
        "original_seconds_per_batch": old,
        "optimized_seconds_per_batch": new,
        "speedup": old / new,
        "repeats": repeats,
        "loops_per_repeat": loops,
    }
    print(
        f"{label:42} {old * 1e3:10.3f} ms -> {new * 1e3:10.3f} ms  {old / new:8.2f}x",
        flush=True,
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_json", type=Path)
    parser.add_argument(
        "--candidate-root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="probe-baseline-") as tmp:
        root = Path(tmp)
        data = json.loads(args.baseline_json.read_text())
        # Extract only the fixed module names needed for these tests.
        for name in (
            "modules/__init__.py",
            "modules/io_tools.py",
            "modules/sequence_utils.py",
            "modules/probe_analysis.py",
        ):
            target = root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(data["files"][name])
        old_su, old_pa = load_modules(root)
        new_su, new_pa = load_modules(args.candidate_root.resolve())

        rng = random.Random(20260919)
        dna = lambda n: "".join(rng.choices("ACGT", k=n))
        probe = dna(25)
        forward, reverse = dna(24), dna(24)
        products = [
            forward + dna(152) + old_su.reverse_complement(reverse) for _ in range(100)
        ]
        exact, three_mm, random_only, iupac = [], [], [], []
        degenerate = list(probe)
        for i in (3, 10, 20):
            degenerate[i] = "N"
        degenerate_probe = "".join(degenerate)
        for index, product in enumerate(products):
            pos = 60 + index % 51
            exact.append(product[:pos] + probe + product[pos + 25 :])
            mutated = list(probe)
            for i in (2, 10, 18):
                mutated[i] = rng.choice([b for b in "ACGT" if b != mutated[i]])
            three_mm.append(product[:pos] + "".join(mutated) + product[pos + 25 :])
            random_only.append(product)
            # Force an ambiguous-target mismatch in a subset of cases.
            iupac.append(
                product[:pos] + ("N" if index % 5 == 0 else "A") + product[pos + 1 :]
            )

        # Directly compare against the uploaded implementation, beyond the
        # independent reference in the unittest suite.
        for index in range(5000):
            p = "".join(
                rng.choices("ACGTRYSWKMBDHVNacgtryswkmbdhvn?!-", k=rng.randrange(0, 36))
            )
            s = "".join(rng.choices("ACGTNRYacgtnry?!-", k=rng.randrange(0, 151)))
            assert old_pa.match_probe(p, s) == new_pa.match_probe(p, s), index
            assert old_su.count_mismatches(p, s) == new_su.count_mismatches(p, s), index
        print(
            "5000 direct randomized comparisons with the uploaded baseline passed.",
            flush=True,
        )

        results = []
        for name, query, dataset in [
            ("Probe search: exact match", probe, exact),
            ("Probe search: three planted mismatches", probe, three_mm),
            ("Probe search: unrelated sequence", probe, random_only),
            ("Probe search: IUPAC query", degenerate_probe, iupac),
        ]:
            results.append(
                time_pair(
                    name,
                    lambda q=query, ds=dataset: [old_pa.match_probe(q, s) for s in ds],
                    lambda q=query, ds=dataset: [new_pa.match_probe(q, s) for s in ds],
                    len(dataset),
                )
            )

        results.append(
            time_pair(
                "Reverse complement: 200 bp",
                lambda: [old_su.reverse_complement(s) for s in products],
                lambda: [new_su.reverse_complement(s) for s in products],
                100,
                loops=50,
            )
        )
        windows = [s[50:75] for s in products]
        results.append(
            time_pair(
                "Mismatch count: 25 bp",
                lambda: [old_su.count_mismatches(probe, s) for s in windows],
                lambda: [new_su.count_mismatches(probe, s) for s in windows],
                100,
                loops=50,
            )
        )

        # Same primers for every product; 30 exact, 50 three-mismatch and
        # 20 unrelated probe regions. All product lengths are 200 bp.
        mixed = exact[:30] + three_mm[30:80] + random_only[80:]
        for duplicate in (False, True):
            records = []
            for i, sequence in enumerate(mixed):
                start = i * 1000 + 100
                coords = f" start={start} end={start + 199} len=200"
                records.append(
                    {"header": f"contig_{2*i+1}" + coords, "sequence": sequence}
                )
                if duplicate:
                    records.append(
                        {
                            "header": f"contig_{2*i+2}" + coords,
                            "sequence": old_su.reverse_complement(sequence),
                        }
                    )
            results.append(
                time_pair(
                    "Normalize products: "
                    + ("paired orientations" if duplicate else "no duplicates"),
                    lambda: old_pa.normalize_ipcr_products(
                        records, forward, reverse, probe, 3
                    ),
                    lambda: new_pa.normalize_ipcr_products(
                        records, forward, reverse, probe, 3
                    ),
                    len(records),
                )
            )

        if hasattr(old_pa, "_parse_ipcr_jsonl_records"):
            for duplicate in (False, True):
                records = []
                for i, sequence in enumerate(mixed):
                    start = i * 1000 + 100
                    record = {
                        "experiment_id": "manual",
                        "sequence_id": "contig-A",
                        "source_file": "synthetic.fna",
                        "start": start,
                        "end": start + 200,
                        "length": 200,
                        "type": "forward",
                        "sequence": sequence,
                    }
                    records.append(record)
                    if duplicate:
                        records.append(
                            dict(
                                record,
                                type="revcomp",
                                sequence=old_su.reverse_complement(sequence),
                            )
                        )
                results.append(
                    time_pair(
                        "Normalize JSONL: "
                        + ("paired orientations" if duplicate else "no duplicates"),
                        lambda: old_pa.normalize_ipcr_products(
                            records, forward, reverse, probe, 3
                        ),
                        lambda: new_pa.normalize_ipcr_products(
                            records, forward, reverse, probe, 3
                        ),
                        len(records),
                    )
                )

        report = {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "seed": 20260919,
            "amplicon_length": 200,
            "probe_length": 25,
            "randomized_direct_equivalence_cases": 5000,
            "notes": "Synthetic in-process timings; median of 7 interleaved repeats; warm query caches; no native engine, filesystem scan, multiprocessing, downloads, or output serialization timed.",
            "results": results,
        }
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
