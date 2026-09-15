from __future__ import annotations

import unittest

from modules.probe_analysis import (
    _parse_fasta_records_from_text,
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

    def test_ipcr_products_are_oriented_and_deduplicated_by_locus(self) -> None:
        forward = "AACCGG"
        reverse = "GGTACC"
        canonical = forward + "TTTTTTTT" + reverse_complement(reverse)
        records = [
            {
                "header": f"manual_1 start=20 end=39 len={len(canonical)}",
                "sequence": canonical,
            },
            {
                "header": f"manual_2 start=20 end=39 len={len(canonical)}",
                "sequence": reverse_complement(canonical),
            },
            {
                "header": f"manual_3 start=200 end=219 len={len(canonical)}",
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
        self.assertEqual(products[0]["product_source_id"], "manual")
        self.assertEqual(products[1]["product_start"], 200)

    def test_unreconcilable_product_fails_loudly(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "cannot be reconciled"):
            normalize_ipcr_products(
                [{"header": "manual_1 start=1 end=20 len=20", "sequence": "A" * 20}],
                "CCCCCC",
                "GGGGGG",
                probe=None,
                max_mismatches=2,
            )


if __name__ == "__main__":
    unittest.main()
