#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diamond_contamination.py

Author: Ankush Sharma (Ankush.Sharma@uga.edu)
Plant Genome Mapping Laboratory, UGA

Description:
  This script screens multiple genomes for contamination using DIAMOND "blastx"
  against a known protein FASTA reference (e.g., UniVec or PhiX proteins).

It reads a TSV file (genome_list.tsv) with two columns:
  GenomeID   /path/to/genome.fna
ignoring lines beginning with '#'.

For each genome in the list:
  - Runs diamond blastx --db <ref>.dmnd --query <genome>.fna
  - Merges hits into one output file, prefixing each line with:
      <GenomeID> <ContaminantLabel> qseqid sseqid pident ...

Usage:
  python diamond_contamination.py \
    --genome_list genome_list.tsv \
    --ref_faa reference_proteins.faa \
    --out univec_compiled.tsv \
    --contaminant_label univec \
    --threads 8

Steps:
  1) Optionally download <ref_faa> if not present (via --download_url).
  2) Makes a diamond DB from <ref_faa>.
  3) Loops over the genome list, runs diamond blastx, compiles all hits.

"""

import sys
import os
import subprocess
import argparse
import logging
from datetime import datetime

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

def maybe_download_ref_faa(ref_faa, url=None):
    """
    If ref_faa doesn't exist and a download URL is given, attempt to wget it.
    """
    if os.path.isfile(ref_faa):
        logging.info(f"Reference FASTA '{ref_faa}' already present.")
        return
    if not url:
        logging.info("No download_url provided; skipping download.")
        return
    logging.info(f"Downloading reference FASTA from {url}")
    cmd = ["wget", "-O", ref_faa, url]
    subprocess.run(cmd, check=True)
    if not os.path.isfile(ref_faa):
        logging.error("Download failed or file missing after download.")
        sys.exit(1)

def make_diamond_db(ref_faa, db_prefix):
    """
    Create a DIAMOND database from the protein FASTA.
    """
    cmd = ["diamond", "makedb", "--in", ref_faa, "-d", db_prefix]
    logging.info("Creating DIAMOND DB:\n  " + " ".join(cmd))
    subprocess.run(cmd, check=True)

def run_diamond_blastx(genome_fasta, db_prefix, threads=1):
    """
    Runs diamond blastx on a single genome, returning lines in outfmt 6.
    """
    cmd = [
        "diamond", "blastx",
        "--db", f"{db_prefix}.dmnd",
        "--query", genome_fasta,
        "--outfmt", "6",
        "--threads", str(threads),
        "--evalue", "1e-5"
    ]
    logging.info(f"Running DIAMOND on {genome_fasta}:\n  {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    lines = result.stdout.strip().split("\n")
    return [l for l in lines if l.strip()]

def main():
    setup_logging()
    parser = argparse.ArgumentParser(description="DIAMOND-based contamination screen for multiple genomes.")
    parser.add_argument("--genome_list", required=True,
                        help="TSV with 2 columns: GenomeID /path/to/genome.fna")
    parser.add_argument("--ref_faa", required=True,
                        help="Protein FASTA of contaminant references (e.g., univec_proteins.faa, phiX_proteins.faa).")
    parser.add_argument("--out", default="compiled_hits.tsv",
                        help="Output TSV with merged hits.")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of CPU threads for DIAMOND.")
    parser.add_argument("--contaminant_label", default="unknown",
                        help="Label added to each line to indicate the reference (e.g. 'univec' or 'phix').")
    parser.add_argument("--download_url", default=None,
                        help="Optional: if ref_faa is missing, download from this URL.")
    args = parser.parse_args()

    start_time = datetime.now()

    # Possibly download ref_faa
    maybe_download_ref_faa(args.ref_faa, args.download_url)

    # Build diamond db
    db_prefix = "tmp_contaminant_db"
    make_diamond_db(args.ref_faa, db_prefix)

    # Validate genome list
    if not os.path.isfile(args.genome_list):
        logging.error(f"Genome list '{args.genome_list}' not found.")
        sys.exit(1)

    # Prepare output, add a header if you wish
    header = "GenomeID\tContaminant\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
    with open(args.out, "w") as fout:
        fout.write(header + "\n")

    total_hits = 0

    # Loop over each line in genome_list
    with open(args.genome_list, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                logging.warning(f"Skipping line (not enough cols): {line}")
                continue

            genome_id = parts[0]
            genome_path = parts[1]
            if not os.path.isfile(genome_path):
                logging.warning(f"Genome file {genome_path} not found for {genome_id}. Skipping.")
                continue

            logging.info(f"=== Processing {genome_id} with reference '{args.contaminant_label}' ===")
            hits = run_diamond_blastx(genome_path, db_prefix, threads=args.threads)
            logging.info(f"Found {len(hits)} hits for {genome_id}")

            total_hits += len(hits)
            with open(args.out, "a") as fout:
                for h in hits:
                    fout.write(f"{genome_id}\t{args.contaminant_label}\t{h}\n")

    elapsed = datetime.now() - start_time
    logging.info(f"Done. Wrote {total_hits} total hits to {args.out}. Elapsed time: {elapsed}.")


if __name__ == "__main__":
    main()
