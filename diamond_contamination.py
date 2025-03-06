#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diamond_contamination.py

Author: Ankush Sharma (Ankush.Sharma@uga.edu)
Plant Genome Mapping Laboratory, UGA

Description:
  This script screens multiple genomes for potential contamination
  using DIAMOND "blastx" searches against a known contaminant protein set.
  (Example: checking for PhiX phage contamination.)

Steps:
  1) Optionally download the contaminant protein FASTA if not present.
  2) Create a DIAMOND DB from that protein FASTA.
  3) Read a simple TSV (genome_list.tsv) with two columns:
       GenomeID  /path/to/genome.fna
  4) For each genome, run "diamond blastx" (translating the genome),
     capture outfmt6 lines, and prefix them with the GenomeID.
  5) Merge all hits into a single output file (default: compiled_contamination.tsv).

Usage (example):
  python diamond_contamination.py \
      --genome_list genome_list.tsv \
      --ref_faa phix_proteins.faa \
      --out compiled_contamination.tsv \
      --threads 8

Requirements:
  - DIAMOND installed and on your $PATH or a loaded module
  - Python 3
"""

import sys
import os
import subprocess
import argparse
import logging
from datetime import datetime

def setup_logging():
    """
    Basic logger for console output with timestamps.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

def maybe_download_ref_faa(ref_faa, url=None):
    """
    If ref_faa is missing and a URL is provided, attempt to download it.
    Adjust or remove if you already have the file.
    """
    if os.path.isfile(ref_faa):
        logging.info(f"{ref_faa} already present. No download needed.")
        return
    if url is None:
        logging.info("No URL provided for reference; not downloading anything.")
        return
    logging.info(f"Downloading reference FASTA from {url}")
    cmd = ["wget", "-O", ref_faa, url]
    subprocess.run(cmd, check=True)
    if not os.path.isfile(ref_faa):
        logging.error(f"Download failed or {ref_faa} not found after download.")
        sys.exit(1)

def make_diamond_db(ref_faa, db_prefix):
    """
    Create a DIAMOND database from the contaminant protein FASTA.
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
    Runs diamond blastx on a single genome.
    Returns the output lines in outfmt6 as a list of strings.
    """
    cmd = [
        "diamond", "blastx",
        "--db", f"{db_prefix}.dmnd",
        "--query", genome_fasta,
        "--outfmt", "6",
        "--threads", str(threads),
        "--evalue", "1e-5"
    ]
    logging.info(f"Running DIAMOND blastx on {genome_fasta}:\n  {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    lines = result.stdout.strip().split("\n")
    # Filter out any empty lines
    return [l for l in lines if l.strip()]

def main():
    setup_logging()

    parser = argparse.ArgumentParser(
        description="Screen multiple genomes for contaminant proteins (via DIAMOND blastx)."
    )
    parser.add_argument("--genome_list", required=True,
                        help="TSV with two columns: GenomeID /path/to/genome.fna")
    parser.add_argument("--ref_faa", required=True,
                        help="Protein FASTA file of the contaminant references (e.g., phiX_proteins.faa).")
    parser.add_argument("--out", default="compiled_contamination.tsv",
                        help="Name of the merged results file.")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of CPU threads to use for diamond blastx.")
    parser.add_argument("--download_url", default=None,
                        help="Optional URL to download ref_faa if not present.")
    args = parser.parse_args()

    start_time = datetime.now()

    # Possibly download the reference if missing
    maybe_download_ref_faa(args.ref_faa, url=args.download_url)

    # Make a diamond DB
    db_prefix = "contaminant_db"
    make_diamond_db(args.ref_faa, db_prefix)

    # Check the genome_list file
    if not os.path.isfile(args.genome_list):
        logging.error(f"Genome list file not found: {args.genome_list}")
        sys.exit(1)

    # Prepare final output
    out_file = args.out
    # Optional header line in compiled file:
    header = "GenomeID\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
    with open(out_file, "w") as fout:
        fout.write(header + "\n")

    total_hits = 0

    # Loop over each genome in the TSV
    with open(args.genome_list, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 2:
                logging.warning(f"Skipping line (not enough columns): {line}")
                continue

            genome_id, genome_path = parts[0], parts[1]
            if not os.path.isfile(genome_path):
                logging.warning(f"Genome file missing for {genome_id}: {genome_path}")
                continue

            logging.info(f"=== Processing {genome_id} ===")
            hits = run_diamond_blastx(genome_path, db_prefix, threads=args.threads)
            logging.info(f"Found {len(hits)} hits for {genome_id}")

            total_hits += len(hits)
            # Write lines to the compiled file, prefixing with GenomeID
            with open(out_file, "a") as fout:
                for h in hits:
                    fout.write(f"{genome_id}\t{h}\n")

    elapsed = datetime.now() - start_time
    logging.info(f"\nAll done. Wrote {total_hits} total hits to {out_file}. Elapsed time: {elapsed}.")

if __name__ == "__main__":
    main()
