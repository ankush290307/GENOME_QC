#!/usr/bin/env python3
"""
Count_pident.py

Author: Ankush Sharma (Ankush.Sharma@uga.edu)
Plant Genome Mapping Laboratory, University of Georgia

Usage:
  python Count_pident.py <blast_outfmt6.tsv> <output_prefix>

Produces:
  1) pident histogram -> <output_prefix>_pident_hist.txt
  2) cumulative distribution -> <output_prefix>_pident_cumulative.txt
  3) prints total queries with hits to stdout
"""

import sys
import math
from collections import defaultdict

def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    blast_file = sys.argv[1]
    out_prefix = sys.argv[2]

    pident_counts = defaultdict(int)
    query_ids = set()

    with open(blast_file, "r") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            qseqid = cols[0]
            try:
                pident = float(cols[2])
            except ValueError:
                continue

            pident_floor = int(math.floor(pident))
            pident_counts[pident_floor] += 1
            query_ids.add(qseqid)

    # Write histogram
    hist_file = f"{out_prefix}_pident_hist.txt"
    with open(hist_file, "w") as out:
        out.write("pident\tcount\n")
        for pid_floor in sorted(pident_counts.keys()):
            out.write(f"{pid_floor}\t{pident_counts[pid_floor]}\n")

    # Write cumulative
    cum_file = f"{out_prefix}_pident_cumulative.txt"
    running_total = 0
    with open(cum_file, "w") as out:
        out.write("pident_threshold\tcumulative_hits\n")
        for pid_floor in sorted(pident_counts.keys(), reverse=True):
            running_total += pident_counts[pid_floor]
            out.write(f">={pid_floor}\t{running_total}\n")

    print(f"\nCount_pident results for: {blast_file}")
    print(f"Total queries with hits: {len(query_ids)}")
    print(f"Histogram: {hist_file}")
    print(f"Cumulative: {cum_file}\n")

if __name__ == "__main__":
    main()

