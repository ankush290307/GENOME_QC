
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diamond_contamination_screen.py

Author: Ankush Sharma (Ankush.Sharma@uga.edu)
Plant Genome Mapping Laboratory, UGA

Description:
  This script screens multiple genomes for possible contaminant sequences
  by using DIAMOND in "blastx" mode. Each genome is translated in 6 frames
  and searched against a known set of contaminant *proteins* (e.g., 
  PhiX protein reference or vector proteins).

Steps:
  1) Reads a two-column TSV file (genome_list.tsv) specifying:
        GenomeID   GenomeFastaPath
     Ignores lines beginning with '#'.
  2) Creates a DIAMOND protein database from your contaminant reference FASTA
     (passed in with --ref_faa).
  3) For each genome in the TSV, runs "diamond blastx --db <contaminant_db> 
     --query <genome> ...", capturing the standard outfmt6 lines.
  4) Prefixes each line of results with the "GenomeID" so you know which
     assembly the hits came from.
  5) Writes a compiled result file (e.g. "compiled_contamination.tsv").

Usage:
  python diamond_contamination_screen.py --genome_list genome_list.tsv \
                                         --ref_faa phix_proteins.faa \
                                         --out compiled_contamination.tsv \
                                         --threads 8

Requirements:
  - DIAMOND >= 0.9
  - Python 3
  - A protein FASTA for your contaminants ( --ref_faa )

Note:
  If you'd rather do a nucleotide-level contamination check (like with 
  UniVec DNA sequences), you typically should use "blastn" instead.

Example:
  1) Prepare phix_proteins.faa or univec_proteins.faa, etc.
  2) Put your genome paths in genome_list.tsv
  3) python diamond_contamination_screen.py \
        --genome_list genome_list.tsv \
        --ref_faa phix_proteins.faa \
        --out phix_contam.tsv \
        --threads 8
"""

import sys
import os
import subprocess
import argparse
import logging
from datetime import datetime

def setup_logging():
    """
    Basic logger that prints to console with timestamps.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

def make_diamond_db(ref_faa, db_prefix):
    """
    Creates a DIAMOND database from the contaminant protein FASTA.
    """
    cmd = [
        "diamond", "makedb",
        "--in", ref_faa,
        "-d", db_prefix
    ]
    logging.info("Creating DIAMOND DB:\n  " + " ".join(cmd))
    subprocess.run(cmd, check=True)

def run_diamond_blastx(genome_fasta, db_prefix, threads=1):
    """
    Runs DIAMOND blastx with outfmt 6, returning the lines of output as a list.
    (No temp file needed; we capture stdout in Python.)
    """
    # Adjust parameters as needed. 
    # Common usage: evalue=1e-5 or 1e-10, etc.
    cmd = [
        "diamond", "blastx",
        "--db", f"{db_prefix}.dmnd",
        "--query", genome_fasta,
        "--outfmt", "6",
        "--threads", str(threads),
        "--evalue", "1e-5"
    ]
    logging.info("Running diamond blastx on " + genome_fasta + ":\n  " + " ".join(cmd))

    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    # result.stdout holds the outfmt6 lines
    lines = result.stdout.strip().split("\n")
    return [x for x in lines if x.strip()]

def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Screen multiple genomes for contamination using diamond blastx."
    )
    parser.add_argument("--genome_list", required=True,
                        help="TSV with two columns: GenomeID, GenomeFastaPath.")
    parser.add_argument("--ref_faa", required=True,
                        help="Protein FASTA of your suspected contaminants (e.g., PhiX, vector proteins).")
    parser.add_argument("--out", default="compiled_contamination.tsv",
                        help="Name of the compiled results TSV.")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of threads for DIAMOND.")
    args = parser.parse_args()

    start_time = datetime.now()

    if not os.path.exists(args.ref_faa):
        logging.error(f"Contaminant reference FASTA not found: {args.ref_faa}")
        sys.exit(1)

    # 1) Create a diamond database from the contaminant proteins
    db_prefix = "contaminant_db"
    make_diamond_db(args.ref_faa, db_prefix)

    # 2) Parse the genome list 
    if not os.path.exists(args.genome_list):
        logging.error("Genome list TSV not found: " + args.genome_list)
        sys.exit(1)

    all_results = []
    with open(args.genome_list, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                logging.warning("Skipping line (not enough columns): " + line)
                continue
            genome_id, genome_path = parts[0], parts[1]
            if not os.path.exists(genome_path):
                logging.warning(f"Genome file does not exist: {genome_path}. Skipping.")
                continue

            logging.info(f"=== Processing genome: {genome_id} ===")
            # 3) Run diamond blastx 
            blast_lines = run_diamond_blastx(genome_path, db_prefix, threads=args.threads)

            # 4) Tag each line with the GenomeID
            # outfmt6 has columns: qseqid sseqid pident length mismatch ...
            # We'll just prefix the entire line with genome_id
            for bl in blast_lines:
                # e.g. "iyGenome1<tab>qseqid sseqid pident length mismatch gapopen..."
                all_results.append(f"{genome_id}\t{bl}")

    # 5) Write compiled file
    logging.info(f"Writing compiled results to {args.out}")
    with open(args.out, "w") as fout:
        # Optionally, you can add a header line 
        # describing the columns. For example:
        # "GenomeID   qseqid  sseqid  pident  length  mismatch  gapopen  qstart..."
        # But you can also just leave it unlabeled.
        for line in all_results:
            fout.write(line + "\n")

    elapsed = datetime.now() - start_time
    logging.info(f"Done. {len(all_results)} total lines. Elapsed time: {elapsed}.")

if __name__ == "__main__":
    main()
