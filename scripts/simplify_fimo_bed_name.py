#!/usr/bin/env python3

import re
import sys
from typing import Optional

# Matches the JASPAR-like motif ID prefix, e.g. "MA0474.4_"
MOTIF_PREFIX_RE = re.compile(r"MA\d+\.\d+_")

def simplify_name(name: str) -> str:
    """
    Input example:
      MER66A_merged_pos_chrY_57199157_57199627_strand_+_MA0474.4_Erg
    Output:
      MER66A|Erg
    """
    # 1) Take subfamily/family before "_merged"
    family = name.split("_merged", 1)[0]

    # 2) Take the TF name after the motif prefix (use the last match if multiple exist)
    matches = list(MOTIF_PREFIX_RE.finditer(name))
    if matches:
        tf = name[matches[-1].end():]
    else:
        # Fallback: if no motif id found, keep everything after last underscore
        tf = name.rsplit("_", 1)[-1]

    return f"{family}|{tf}"

def main(in_bed: str, out_bed: str) -> None:
    with open(in_bed, "r", encoding="utf-8") as fin, open(out_bed, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                fout.write(line)
                continue

            fields[3] = simplify_name(fields[3])
            fout.write("\t".join(fields) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage: simplify_bed_name.py <input.bed> <output.bed>\n"
        )
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
