"""Behavioral checks for the performance-only changes (no timing assertions)."""

from __future__ import annotations

import itertools
import json
import random
import subprocess
from pathlib import Path
import unittest
from unittest.mock import patch

from modules import probe_analysis as pa
from modules import sequence_utils as su

# A deliberately simple reference independent of the optimized lookup table.
_ALLOWED = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


def reference_count(primer: str, window: str) -> int:
    return sum(
        g not in _ALLOWED.get(p, "") for p, g in zip(primer.upper(), window.upper())
    ) + abs(len(primer) - len(window))


def reference_match(probe: str, sequence: str) -> tuple[int, int]:
    if not probe or not sequence or len(sequence) < len(probe):
        return -1, -1
    best, pos = len(probe), -1
    for i in range(len(sequence) - len(probe) + 1):
        mm = reference_count(probe, sequence[i : i + len(probe)])
        if mm < best:
            best, pos = mm, i
            if best == 0:
                break
    return best, pos


class SequencePerformanceRegressions(unittest.TestCase):
    def test_reverse_complement_all_symbols(self):
        before = "ACGTRYSWKMBDHVNacgtryswkmbdhvn-?!ßſ🙂"
        after = "TGCAYRSWMKVHDBNtgcayrswmkvhdbn-?!ßſ🙂"
        self.assertEqual(su.reverse_complement(before), after[::-1])
        self.assertEqual(su.reverse_complement(""), "")

    def test_count_all_iupac_and_genome_symbols(self):
        for p, g in itertools.product(
            "ACGTRYSWKMBDHVNacgtryswkmbdhvn?-ßſ", "ACGTNRYSWKMBDHVNacgtn?-ßſ"
        ):
            self.assertEqual(su.count_mismatches(p, g), reference_count(p, g), (p, g))

    def test_unequal_lengths(self):
        for p, s in [
            ("", "ACG"),
            ("ACGN", ""),
            ("AR", "AGGG"),
            ("ACGT", "AT"),
            ("ß", "GG"),
        ]:
            self.assertEqual(su.count_mismatches(p, s), reference_count(p, s))

    def test_invalid_or_missing_windows(self):
        for p, s in [("", "ACGT"), ("ACGT", ""), ("AAA", "AA")]:
            self.assertEqual(pa.match_probe(p, s), (-1, -1))

    def test_all_mismatch_sentinel_is_preserved(self):
        self.assertEqual(pa.match_probe("AAA", "CCCC"), (3, -1))
        self.assertEqual(pa.match_probe("NN", "NNN"), (2, -1))

    def test_genome_n_is_not_a_wildcard(self):
        self.assertEqual(pa.match_probe("AN", "AN"), (1, 0))
        self.assertEqual(pa.match_probe("N", "NA"), (0, 1))
        self.assertEqual(pa.match_probe("AN", "ATAN"), (0, 0))

    def test_leftmost_ties(self):
        self.assertEqual(pa.match_probe("AC", "TTACAC"), (0, 2))
        self.assertEqual(pa.match_probe("AA", "ACAG"), (1, 0))
        self.assertEqual(pa.match_probe("R", "NCAG"), (0, 2))

    def test_lowercase_and_ambiguous_queries(self):
        for p, s in [
            ("arYn", "TTACTAACGT"),
            ("acgt", "nnACGTacgt"),
            ("NRY", "NACGT"),
            ("WKMBDHV", "ACGTACGTACGT"),
        ]:
            self.assertEqual(pa.match_probe(p, s), reference_match(p, s))

    def test_unicode_retains_original_windowing(self):
        for p, s in [
            ("ß", "GG"),
            ("ſ", "CG"),
            ("aß", "ACGAG"),
            ("AC", "aßac"),
            ("NN", "ß"),
        ]:
            self.assertEqual(pa.match_probe(p, s), reference_match(p, s), (p, s))

    def test_random_equivalence_5000_cases(self):
        rng = random.Random(20260919)
        query_alphabet = "ACGTRYSWKMBDHVNacgtryswkmbdhvn?!-"
        target_alphabet = "ACGTNRYacgtnry?!-"
        for case in range(5000):
            p = "".join(rng.choices(query_alphabet, k=rng.randrange(0, 36)))
            s = "".join(rng.choices(target_alphabet, k=rng.randrange(0, 151)))
            self.assertEqual(su.count_mismatches(p, s), reference_count(p, s), case)
            self.assertEqual(pa.match_probe(p, s), reference_match(p, s), case)

    def test_exhaustive_short_queries(self):
        for probe in map("".join, itertools.product("ARN", repeat=2)):
            for seq in map("".join, itertools.product("ACGN", repeat=4)):
                self.assertEqual(
                    pa.match_probe(probe, seq),
                    reference_match(probe, seq),
                    (probe, seq),
                )

    def test_query_cache_is_bounded(self):
        su._compile_primer.cache_clear()
        for i in range(300):
            su.count_mismatches("A" + str(i), "ACGT")
        self.assertLessEqual(su._compile_primer.cache_info().currsize, 128)


class DeduplicationPerformanceRegressions(unittest.TestCase):
    forward = "AACCGG"
    reverse = "GGTACC"

    def records(self, duplicate=True, coordinate=True):
        product = self.forward + "TTTTTTTT" + su.reverse_complement(self.reverse)
        coords = " start=20 end=39 len=20" if coordinate else ""
        records = [{"header": "manual_1" + coords, "sequence": product}]
        if duplicate:
            records.append(
                {
                    "header": "manual_2" + coords,
                    "sequence": su.reverse_complement(product),
                }
            )
        return records

    def test_duplicate_is_probed_only_once(self):
        with patch.object(pa, "match_probe", wraps=pa.match_probe) as match:
            result = pa.normalize_ipcr_products(
                self.records(), self.forward, self.reverse, "TTTT", 2
            )
        self.assertEqual(match.call_count, 1)
        self.assertEqual(result[0]["duplicate_record_count"], 2)
        self.assertEqual(result[0]["product_id"], "manual_1")
        self.assertEqual(result[0]["product_orientation"], "forward")

    def test_missing_coordinates_keep_both_records(self):
        with patch.object(pa, "match_probe", wraps=pa.match_probe) as match:
            result = pa.normalize_ipcr_products(
                self.records(coordinate=False), self.forward, self.reverse, "TTTT", 2
            )
        self.assertEqual(len(result), 2)
        self.assertEqual(match.call_count, 2)

    def test_identical_sequence_different_loci_keep_both(self):
        records = self.records()
        records[1]["header"] = "manual_2 start=200 end=219 len=20"
        result = pa.normalize_ipcr_products(
            records, self.forward, self.reverse, "TTTT", 2
        )
        self.assertEqual(len(result), 2)

    def test_conflicting_duplicate_still_fails(self):
        records = self.records(duplicate=False)
        other = self.forward + "TTTTTTTA" + su.reverse_complement(self.reverse)
        records.append({"header": "manual_2 start=20 end=39 len=20", "sequence": other})
        with self.assertRaisesRegex(RuntimeError, "conflicting normalized sequences"):
            pa.normalize_ipcr_products(records, self.forward, self.reverse, "TTTT", 2)

    def test_invalid_duplicate_still_fails_primer_validation(self):
        records = self.records(duplicate=False)
        records.append(
            {"header": "manual_2 start=20 end=39 len=20", "sequence": "A" * 20}
        )
        with self.assertRaisesRegex(RuntimeError, "cannot be reconciled"):
            pa.normalize_ipcr_products(records, self.forward, self.reverse, "TTTT", 2)

    def test_no_probe_retains_null_fields(self):
        result = pa.normalize_ipcr_products(
            self.records(), self.forward, self.reverse, None, 2
        )
        self.assertIsNone(result[0]["probe_mismatches"])
        self.assertIsNone(result[0]["probe_position"])


class IPCRExitStatusRegressions(unittest.TestCase):
    def run_status(self, status, stdout=b"", stderr=b""):
        completed = subprocess.CompletedProcess(["ipcr"], status, stdout, stderr)
        with patch.object(pa.io_tools, "run_command", return_value=completed) as run:
            try:
                return pa._run_ipcr_records(
                    "ACGT", "TGCA", Path("genome.fna"), min_len=60, max_len=200
                )
            finally:
                self.assertIs(run.call_args.kwargs["check"], False)

    def test_no_hits_is_empty_not_an_exception(self):
        self.assertEqual(self.run_status(1), [])

    def test_success_still_parses_records(self):
        payload = {
            "experiment_id": "manual",
            "sequence_id": "contig-A",
            "source_file": "genome.fna",
            "start": 1,
            "end": 5,
            "length": 4,
            "type": "forward",
            "seq": "acgt",
        }
        records = self.run_status(0, (json.dumps(payload) + "\n").encode())
        self.assertEqual(
            records,
            [
                {
                    "experiment_id": "manual",
                    "sequence_id": "contig-A",
                    "source_file": "genome.fna",
                    "start": 1,
                    "end": 5,
                    "length": 4,
                    "type": "forward",
                    "sequence": "ACGT",
                    "fwd_mm": 0,
                    "rev_mm": 0,
                }
            ],
        )

    def test_success_without_products_is_empty(self):
        self.assertEqual(self.run_status(0), [])

    def test_jsonl_command_is_preserved(self):
        with patch.object(pa.io_tools, "run_command", return_value=None) as run:
            result = pa._run_ipcr_records(
                "ACGT",
                "TGCA",
                Path("genome.fna"),
                min_len=60,
                max_len=200,
                dry_run=True,
            )
        command = run.call_args.args[0]
        self.assertEqual(command[command.index("--output") + 1], "jsonl")
        self.assertIn("--products", command)
        self.assertIs(run.call_args.kwargs["check"], False)
        self.assertIs(run.call_args.kwargs["dry_run"], True)
        self.assertEqual(result, [])

    def test_malformed_jsonl_still_raises(self):
        with self.assertRaisesRegex(ValueError, "not valid JSON"):
            self.run_status(0, b">not-jsonl\nACGT\n")

    def test_configuration_failure_still_raises(self):
        with self.assertRaisesRegex(ValueError, "usage/config error.*bad arguments"):
            self.run_status(2, stderr=b"bad arguments")

    def test_other_failure_still_raises(self):
        with self.assertRaisesRegex(RuntimeError, "exit 7.*broken engine"):
            self.run_status(7, stderr=b"broken engine")

    def test_missing_executable_still_raises(self):
        with patch.object(
            pa.io_tools, "run_command", side_effect=FileNotFoundError("missing")
        ):
            with self.assertRaises(FileNotFoundError):
                pa._run_ipcr_records(
                    "ACGT", "TGCA", Path("genome.fna"), min_len=60, max_len=200
                )


class JSONLDeduplicationRegressions(unittest.TestCase):
    forward = "AACCGG"
    reverse = "GGTACC"

    def record(self, **changes):
        sequence = self.forward + "TTTTTTTT" + su.reverse_complement(self.reverse)
        record = {
            "experiment_id": "manual",
            "sequence_id": "contig-A",
            "source_file": "genome.fna",
            "start": 20,
            "end": 40,
            "length": len(sequence),
            "type": "forward",
            "sequence": sequence,
        }
        record.update(changes)
        return record

    def normalize(self, records, probe="TTTT"):
        return pa.normalize_ipcr_products(
            records,
            self.forward,
            self.reverse,
            probe,
            max_mismatches=2,
        )

    def test_jsonl_duplicate_is_probed_once_and_first_metadata_is_retained(self):
        first = self.record()
        duplicate = self.record(
            type="revcomp",
            sequence=su.reverse_complement(first["sequence"]),
        )
        with patch.object(pa, "match_probe", wraps=pa.match_probe) as match:
            products = self.normalize([first, duplicate])
        self.assertEqual(match.call_count, 1)
        self.assertEqual(len(products), 1)
        item = products[0]
        self.assertEqual(item["product"], first["sequence"])
        self.assertEqual(item["duplicate_record_count"], 2)
        self.assertEqual(item["product_id"], "manual:contig-A:20-40:forward")
        self.assertEqual(
            item["product_header"],
            "manual:contig-A:20-40:forward source_file=genome.fna len=20",
        )
        self.assertEqual(item["product_source_id"], "contig-A")
        self.assertEqual(item["product_orientation"], "forward")
        self.assertEqual(item["probe_mismatches"], 0)
        self.assertEqual(item["probe_position"], 6)

    def test_jsonl_same_coordinates_on_distinct_contigs_are_not_collapsed(self):
        with patch.object(pa, "match_probe", wraps=pa.match_probe) as match:
            products = self.normalize(
                [self.record(), self.record(sequence_id="contig-B")]
            )
        self.assertEqual(len(products), 2)
        self.assertEqual(match.call_count, 2)
        self.assertEqual([p["duplicate_record_count"] for p in products], [1, 1])

    def test_jsonl_same_coordinates_in_distinct_source_files_are_not_collapsed(self):
        products = self.normalize([self.record(), self.record(source_file="other.fna")])
        self.assertEqual(len(products), 2)

    def test_jsonl_same_sequence_at_distinct_loci_is_not_collapsed(self):
        products = self.normalize([self.record(), self.record(start=120, end=140)])
        self.assertEqual(len(products), 2)

    def test_jsonl_swapped_coordinates_keep_first_metadata(self):
        products = self.normalize([self.record(), self.record(start=40, end=20)])
        self.assertEqual(len(products), 1)
        self.assertEqual(products[0]["duplicate_record_count"], 2)
        self.assertEqual(products[0]["product_start"], 20)
        self.assertEqual(products[0]["product_end"], 40)

    def test_jsonl_conflicting_duplicate_still_fails(self):
        other = self.forward + "TTTTTTTA" + su.reverse_complement(self.reverse)
        with self.assertRaisesRegex(RuntimeError, "conflicting normalized sequences"):
            self.normalize([self.record(), self.record(sequence=other)])

    def test_jsonl_invalid_duplicate_still_fails_primer_validation(self):
        with self.assertRaisesRegex(RuntimeError, "cannot be reconciled"):
            self.normalize([self.record(), self.record(sequence="A" * 20)])

    def test_jsonl_without_probe_retains_null_fields(self):
        with patch.object(pa, "match_probe", wraps=pa.match_probe) as match:
            products = self.normalize([self.record(), self.record()], probe=None)
        match.assert_not_called()
        self.assertIsNone(products[0]["probe_mismatches"])
        self.assertIsNone(products[0]["probe_position"])

    def test_jsonl_and_legacy_locus_keys_remain_distinct(self):
        first = self.record()
        legacy = {
            "header": "contig-A_1 start=20 end=40 len=20",
            "sequence": first["sequence"],
        }
        products = self.normalize([first, legacy])
        self.assertEqual(len(products), 2)


if __name__ == "__main__":
    unittest.main()
