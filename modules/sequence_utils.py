"""
modules/sequence_utils.py

Basic utilities for sequence analysis.
"""

def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence, handling ambiguous bases and case.
    """
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
        'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h',
        'h': 'd', 'v': 'b', 'n': 'n'
    }
    return ''.join(complement.get(base, base) for base in reversed(seq))

# Bit masks: A=1, C=2, G=4, T=8
_IUPAC = {
    "A": 0b0001, "C": 0b0010, "G": 0b0100, "T": 0b1000,
    "R": 0b0101, "Y": 0b1010, "S": 0b0110, "W": 0b1001,
    "K": 0b1100, "M": 0b0011, "B": 0b1110, "D": 0b1101,
    "H": 0b1011, "V": 0b0111, "N": 0b1111,
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

def count_mismatches(primer: str, window: str) -> int:
    """
    Count mismatches between primer (may contain IUPAC codes) and a genome window.
    Args:
      primer: primer/probe sequence (IUPAC allowed)
      window: genome subsequence (expected A/C/G/T/N)
    """
    mm = 0
    for p, g in zip(primer.upper(), window.upper()):
        if not _base_match(g, p):  # NOTE: genome first, primer second (fixed)
            mm += 1
    mm += abs(len(primer) - len(window))
    return mm


# ---