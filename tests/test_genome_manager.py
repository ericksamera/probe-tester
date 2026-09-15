from __future__ import annotations

import json
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from modules.genome_manager import get_genomes_mapping


class GenomeManagerTests(unittest.TestCase):
    def test_genbank_mapping_does_not_require_asm_assembly_name(self) -> None:
        payload = {
            "reports": [
                {
                    "current_accession": "GCA_000089865.1",
                    "organism": {"organism_name": "Mycoplasmopsis agalactiae"},
                    "assembly_info": {"assembly_name": "ASM8986v1"},
                },
                {
                    "current_accession": "GCA_900088695.1",
                    "organism": {"organism_name": "Mycoplasmopsis agalactiae"},
                    "assembly_info": {"assembly_name": "JF4428"},
                },
                {
                    "current_accession": "GCF_900088695.1",
                    "organism": {"organism_name": "Mycoplasmopsis agalactiae"},
                    "assembly_info": {"assembly_name": "ASM90008869v1"},
                },
                {
                    "current_accession": None,
                    "organism": {"organism_name": "Mycoplasmopsis agalactiae"},
                    "assembly_info": {"assembly_name": "ASM_missing"},
                },
            ]
        }
        completed = SimpleNamespace(stdout=json.dumps(payload).encode())

        with patch("modules.genome_manager.io_tools.run_command", return_value=completed) as run:
            mapping = get_genomes_mapping(2110)

        run.assert_called_once_with(
            [
                "datasets",
                "summary",
                "genome",
                "taxon",
                "2110",
                "--assembly-source",
                "GenBank",
                "--assembly-version",
                "latest",
            ],
            capture_output=True,
        )
        self.assertEqual(
            mapping,
            {
                "Mycoplasmopsis-agalactiae": [
                    "GCA_000089865.1",
                    "GCA_900088695.1",
                ]
            },
        )


if __name__ == "__main__":
    unittest.main()
