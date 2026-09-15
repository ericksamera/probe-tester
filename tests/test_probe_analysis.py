from __future__ import annotations

import json
import unittest

from modules.probe_analysis import (
    _parse_fasta_records_from_text,
    _parse_ipcr_jsonl_records,
    normalize_ipcr_products,
)
from modules.sequence_utils import reverse_complement


class ProbeAnalysisTests(unittest.TestCase):
    def test_fasta_parser_preserves_headers(self) -> None:
        records = _parse_fasta_records_from_text(
            ">manual_1 start=10 end=89 len=80\nACGT\n"
            ">manual_2 start=100 end=179 len=80\nTGCA\n"
        )
        self.assertEqual(records[0]["header"], "manual_1 start=10 end=89 len=80")
        self.assertEqual(records[0]["sequence"], "ACGT")
        self.assertEqual(records[1]["header"], "manual_2 start=100 end=179 len=80")

    def test_ipcr_jsonl_parser_preserves_sequence_identity(self) -> None:
        payload = {
            "experiment_id": "manual",
            "sequence_id": "contig-A",
            "start": 20,
            "end": 40,
            "length": 20,
            "type": "forward",
            "seq": "ACGT" * 5,
            "source_file": "genome.fna",
        }
        records = _parse_ipcr_jsonl_records(json.dumps(payload) + "\n")
        self.assertEqual(records[0]["sequence_id"], "contig-A")
        self.assertEqual(records[0]["source_file"], "genome.fna")
        self.assertEqual(records[0]["sequence"], "ACGT" * 5)
        self.assertEqual(records[0]["fwd_mm"], 0)
        self.assertEqual(records[0]["rev_mm"], 0)

    def test_ipcr_products_are_oriented_and_deduplicated_by_locus(self) -> None:
        forward = "AACCGG"
        reverse = "GGTACC"
        canonical = forward + "TTTTTTTT" + reverse_complement(reverse)
        records = [
            {
                "experiment_id": "manual",
                "sequence_id": "contig-A",
                "source_file": "genome.fna",
                "start": 20,
                "end": 40,
                "length": len(canonical),
                "type": "forward",
                "sequence": canonical,
            },
            {
                "experiment_id": "manual",
                "sequence_id": "contig-A",
                "source_file": "genome.fna",
                "start": 20,
                "end": 40,
                "length": len(canonical),
                "type": "revcomp",
                "sequence": reverse_complement(canonical),
            },
            {
                "experiment_id": "manual",
                "sequence_id": "contig-B",
                "source_file": "genome.fna",
                "start": 20,
                "end": 40,
                "length": len(canonical),
                "type": "forward",
                "sequence": canonical,
            },
        ]
        products = normalize_ipcr_products(
            records,
            forward,
            reverse,
            probe="TTTT",
            max_mismatches=2,
        )
        self.assertEqual(len(products), 2)
        self.assertEqual(products[0]["product"], canonical)
        self.assertEqual(products[0]["forward_mismatches"], 0)
        self.assertEqual(products[0]["reverse_mismatches"], 0)
        self.assertEqual(products[0]["duplicate_record_count"], 2)
        self.assertEqual(products[0]["product_source_id"], "contig-A")
        self.assertEqual(products[1]["product_source_id"], "contig-B")

    def test_legacy_fasta_records_still_normalize(self) -> None:
        forward = "AACCGG"
        reverse = "GGTACC"
        canonical = forward + "TTTTTTTT" + reverse_complement(reverse)
        products = normalize_ipcr_products(
            [
                {
                    "header": f"manual_1 start=20 end=40 len={len(canonical)}",
                    "sequence": canonical,
                }
            ],
            forward,
            reverse,
            probe="TTTT",
            max_mismatches=2,
        )
        self.assertEqual(len(products), 1)
        self.assertEqual(products[0]["product_source_id"], "manual")

    def test_unreconcilable_product_fails_loudly(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "cannot be reconciled"):
            normalize_ipcr_products(
                [
                    {
                        "experiment_id": "manual",
                        "sequence_id": "contig-A",
                        "source_file": "genome.fna",
                        "start": 1,
                        "end": 21,
                        "length": 20,
                        "type": "forward",
                        "sequence": "A" * 20,
                    }
                ],
                "CCCCCC",
                "GGGGGG",
                probe=None,
                max_mismatches=2,
            )


if __name__ == "__main__":
    unittest.main()
