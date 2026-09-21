from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from modules import marker_manager as mm
from modules.fasta_io import read_fasta


class MarkerManagerTests(unittest.TestCase):
    def test_download_bold_marker_groups_species_and_filters_marker(self) -> None:
        preprocess_payload = json.dumps(
            {
                "successful_terms": [
                    {
                        "submitted": "tax:Ixodes",
                        "matched": "tax:genus:Ixodes,12345",
                    }
                ]
            }
        )
        query_payload = json.dumps({"query_id": "abc=="})
        tsv = (
            "processid\tmarker_code\tspecies\tidentification\tnuc\tinsdc_acs\tbin_uri\n"
            "P1\tCOI-5P\tIxodes scapularis\tIxodes scapularis\tAC-GT N\tK1\tBOLD:A\n"
            "P2\t12S\tIxodes scapularis\tIxodes scapularis\tACGT\tK2\t\n"
            "P3\tCOI-5P\tIxodes pacificus\tIxodes pacificus\tAACCGG\t\t\n"
            "P4\tCOI-5P\t\tIxodes sp.\tAACCGG\t\t\n"
        )

        with tempfile.TemporaryDirectory() as tmp, patch.object(
            mm, "_get_text", side_effect=[preprocess_payload, query_payload, tsv]
        ):
            manifest = mm.download_bold_marker(
                "Ixodes", marker="COI", outdir=Path(tmp), min_length=4
            )
            base = Path(tmp) / "markers" / "COI-5P"
            scap = list(read_fasta(base / "Ixodes-scapularis" / "sequences.fasta"))
            pac = list(read_fasta(base / "Ixodes-pacificus" / "sequences.fasta"))

        self.assertEqual(manifest["kept_records"], 2)
        self.assertEqual(manifest["skipped_marker"], 1)
        self.assertEqual(manifest["skipped_unidentified"], 1)
        self.assertEqual(scap[0][0], "P1")
        self.assertEqual(scap[0][2], "ACGTN")
        self.assertEqual(pac[0][0], "P3")

    def test_bold_query_preprocesses_taxon_before_query(self) -> None:
        preprocess_payload = {
            "successful_terms": [
                {
                    "submitted": "tax:Penaeidae",
                    "matched": "tax:family:Penaeidae,7127",
                }
            ]
        }
        with patch.object(
            mm, "_get_json", side_effect=[preprocess_payload, {"query_id": "qid=="}]
        ) as get_json:
            query_id = mm._bold_query_id("Penaeidae")

        self.assertEqual(query_id, "qid==")
        first_url = get_json.call_args_list[0].args[0]
        second_url = get_json.call_args_list[1].args[0]
        self.assertIn("/query/preprocessor?", first_url)
        self.assertIn("query=tax%3APenaeidae", first_url)
        self.assertIn("/query?", second_url)
        self.assertIn("query=tax%3Afamily%3APenaeidae", second_url)
        self.assertIn("extent=full", second_url)

    def test_bold_query_rejects_unresolved_taxon(self) -> None:
        with patch.object(mm, "_get_json", return_value={"successful_terms": []}):
            with self.assertRaisesRegex(ValueError, "could not resolve taxon"):
                mm._bold_query_id("not-a-real-taxon")

    def test_marker_aliases(self) -> None:
        self.assertEqual(mm.normalize_marker_name("COI"), "COI-5P")
        self.assertEqual(mm.normalize_marker_name("co1"), "COI-5P")
        self.assertEqual(mm.normalize_marker_name("12S"), "12S")


if __name__ == "__main__":
    unittest.main()
