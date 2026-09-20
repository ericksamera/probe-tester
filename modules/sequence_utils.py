"""
modules/sequence_utils.py

Basic utilities for sequence analysis.
"""

from functools import lru_cache
from typing import Tuple

_COMPLEMENT_TABLE = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


def reverse_complement(seq: str) -> str:
    """Reverse-complement DNA, retaining case and passing unknown symbols through."""
    return seq.translate(_COMPLEMENT_TABLE)[::-1]


# Bit masks: A=1, C=2, G=4, T=8
_IUPAC = {
    "A": 0b0001,
    "C": 0b0010,
    "G": 0b0100,
    "T": 0b1000,
    "R": 0b0101,
    "Y": 0b1010,
    "S": 0b0110,
    "W": 0b1001,
    "K": 0b1100,
    "M": 0b0011,
    "B": 0b1110,
    "D": 0b1101,
    "H": 0b1011,
    "V": 0b0111,
    "N": 0b1111,
}


def _base_match(genome_base: str, primer_base: str) -> bool:
    """
    True if primer_base (IUPAC) can match genome_base (A/C/G/T only).
    Genome 'N' (or non-ACGT) is treated as a mismatch, mirroring ipcr. :contentReference[oaicite:0]{index=0}
    """
    g = genome_base.upper()
    p = primer_base.upper()
    if g not in "ACGT":
        return False
    return (_IUPAC.get(p, 0) & _IUPAC[g]) != 0


# Only unambiguous genomic bases can match. In particular, genome N is
# still a mismatch even when the primer/probe contains N.
_GENOME_MASKS = {"A": 0b0001, "C": 0b0010, "G": 0b0100, "T": 0b1000}


@lru_cache(maxsize=128)
def _compile_primer(primer_upper: str) -> Tuple[Tuple[int, ...], bool]:
    """Cache query encodings in a bounded cache, never target sequences."""
    masks = tuple(_IUPAC.get(base, 0) for base in primer_upper)
    unambiguous = all(mask in (1, 2, 4, 8) for mask in masks)
    return masks, unambiguous


def count_mismatches(primer: str, window: str) -> int:
    """Count IUPAC-aware mismatches, preserving unequal-length behavior."""
    masks, _ = _compile_primer(primer.upper())
    genome_mask = _GENOME_MASKS.get
    return sum(
        not (primer_mask & genome_mask(base, 0))
        for primer_mask, base in zip(masks, window.upper())
    ) + abs(len(primer) - len(window))
