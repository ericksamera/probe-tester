"""Download taxonomically scoped marker references from BOLD."""

from __future__ import annotations

import csv
import io
import json
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

from modules.fasta_io import write_fasta

logger = logging.getLogger(__name__)

BOLD_API_BASE = "https://portal.boldsystems.org/api"
_BOLD_FIELDS = (
    "processid",
    "marker_code",
    "species",
    "identification",
    "nuc",
    "insdc_acs",
    "bin_uri",
)
_MARKER_ALIASES = {
    "COI": "COI-5P",
    "CO1": "COI-5P",
    "COI-5P": "COI-5P",
}


def normalize_marker_name(marker: str) -> str:
    value = marker.strip().upper()
    return _MARKER_ALIASES.get(value, marker.strip())


def _get_text(url: str, *, timeout: int = 120) -> str:
    request = Request(
        url,
        headers={
            "Accept": "application/json, text/tab-separated-values;q=0.9, */*;q=0.1",
            "User-Agent": "probe-tester/marker-downloader",
        },
    )
    with urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8-sig")


def _get_json(url: str, *, timeout: int = 120) -> dict:
    payload = json.loads(_get_text(url, timeout=timeout))
    if not isinstance(payload, dict):
        raise ValueError("BOLD API returned a non-object JSON response")
    return payload


def _bold_query_id(taxon: str, *, api_base: str = BOLD_API_BASE) -> str:
    params = urlencode({"query": f"tax:{taxon}", "extent": "full"})
    payload = _get_json(f"{api_base.rstrip('/')}/query?{params}")
    query_id = payload.get("query_id")
    if not isinstance(query_id, str) or not query_id:
        raise ValueError("BOLD API query response did not contain query_id")
    return query_id


def _download_bold_tsv(query_id: str, *, api_base: str = BOLD_API_BASE) -> str:
    fields = ",".join(_BOLD_FIELDS)
    encoded_id = quote(query_id, safe="")
    params = urlencode({"format": "tsv", "fields": fields})
    return _get_text(
        f"{api_base.rstrip('/')}/documents/{encoded_id}/download?{params}",
        timeout=300,
    )


def _clean_sequence(sequence: str) -> str:
    # BOLD exports can contain alignment padding. ipcr expects raw sequence.
    return re.sub(r"[^A-Za-z]", "", sequence or "").upper()


def _species_slug(species: str) -> str:
    return (
        re.sub(r"[^A-Za-z0-9._-]+", "-", species.strip()).strip("-") or "unidentified"
    )


def _iter_bold_records(tsv_text: str) -> Iterable[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(tsv_text), delimiter="\t")
    if not reader.fieldnames:
        return
    for row in reader:
        yield {str(k): (v or "") for k, v in row.items() if k is not None}


def download_bold_marker(
    taxon: str,
    *,
    marker: str = "COI-5P",
    outdir: Path = Path("."),
    min_length: int = 0,
    max_records: Optional[int] = None,
    api_base: str = BOLD_API_BASE,
) -> dict:
    """Download one BOLD marker for a taxon into species-level multi-record FASTAs."""
    marker = normalize_marker_name(marker)
    query_id = _bold_query_id(taxon, api_base=api_base)
    tsv_text = _download_bold_tsv(query_id, api_base=api_base)

    by_species: Dict[str, List[tuple[str, str, str]]] = defaultdict(list)
    seen_ids: set[str] = set()
    stats = {
        "downloaded_rows": 0,
        "kept_records": 0,
        "skipped_marker": 0,
        "skipped_unidentified": 0,
        "skipped_short": 0,
        "skipped_missing_sequence": 0,
        "skipped_duplicate_id": 0,
    }

    for row in _iter_bold_records(tsv_text):
        stats["downloaded_rows"] += 1
        if row.get("marker_code", "").strip().upper() != marker.upper():
            stats["skipped_marker"] += 1
            continue

        species = row.get("species", "").strip()
        if not species:
            stats["skipped_unidentified"] += 1
            continue

        sequence = _clean_sequence(row.get("nuc", ""))
        if not sequence:
            stats["skipped_missing_sequence"] += 1
            continue
        if len(sequence) < min_length:
            stats["skipped_short"] += 1
            continue

        processid = row.get("processid", "").strip()
        record_id = processid or row.get("insdc_acs", "").strip()
        if not record_id:
            record_id = f"bold_record_{stats['downloaded_rows']}"
        if record_id in seen_ids:
            stats["skipped_duplicate_id"] += 1
            continue
        seen_ids.add(record_id)

        metadata = [
            f"species={species.replace(' ', '_')}",
            f"marker={marker}",
        ]
        if row.get("insdc_acs"):
            metadata.append(f"insdc={row['insdc_acs'].strip()}")
        if row.get("bin_uri"):
            metadata.append(f"bin={row['bin_uri'].strip()}")

        by_species[species].append((record_id, " ".join(metadata), sequence))
        stats["kept_records"] += 1
        if max_records is not None and stats["kept_records"] >= max_records:
            break

    marker_dir = Path(outdir) / "markers" / marker
    marker_dir.mkdir(parents=True, exist_ok=True)

    for species, records in sorted(by_species.items()):
        species_dir = marker_dir / _species_slug(species)
        species_dir.mkdir(parents=True, exist_ok=True)
        write_fasta(records, species_dir / "sequences.fasta")

    manifest = {
        "source": "BOLD",
        "taxon": taxon,
        "marker": marker,
        "query_id": query_id,
        "species_count": len(by_species),
        **stats,
    }
    (marker_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    logger.info(
        "Downloaded %d %s records across %d species into %s",
        stats["kept_records"],
        marker,
        len(by_species),
        marker_dir,
    )
    return manifest
