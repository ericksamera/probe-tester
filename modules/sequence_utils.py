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

def count_mismatches(a: str, b: str) -> int:
    """
    Count the number of mismatches between two sequences of possibly different lengths.
    """
    return sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))

# ---