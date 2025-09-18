"""
modules/formatting.py
"""
from typing import List

def format_table(headers: List[str], rows: List[str], style="plain", padding=1) -> str:
    """
    Pretty-print a table with aligned columns (no third-party deps).
    - style: "grid" (ASCII box) or "plain" (no borders)
    - numbers are right-aligned; text is left-aligned
    """
    import re
    H = [str(h) for h in headers]
    R = [[("" if c is None else str(c)) for c in row] for row in rows]
    ncols = len(H)
    num_re = re.compile(r"^[+-]?(?:\d+(?:\.\d+)?|\.\d+)%?$")

    # column widths
    widths = [len(H[i]) for i in range(ncols)]
    for row in R:
        for i, cell in enumerate(row):
            if i < ncols:
                widths[i] = max(widths[i], len(cell))

    def is_num(s): return bool(num_re.match(s))
    align = ['>' if all(is_num(r[i]) for r in R if i < len(r)) else '<' for i in range(ncols)]

    pad = " " * padding
    def fmt_row(cells):
        cells = [(cells[i] if i < len(cells) else "") for i in range(ncols)]
        return "|" + "|".join(f"{pad}{cells[i]:{align[i]}{widths[i]}}{pad}" for i in range(ncols)) + "|"

    if style == "grid":
        sep = "+" + "+".join("-" * (w + 2 * padding) for w in widths) + "+"
        out = [sep, fmt_row(H), sep]
        out += [fmt_row(r) for r in R]
        out.append(sep)
        return "\n".join(out)
    else:  # plain
        out = [fmt_row(H)]
        out += [fmt_row(r) for r in R]
        return "\n".join(out)
