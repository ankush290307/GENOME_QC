#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
busco_diamond.py

Author: Ankush Sharma (Ankush.Sharma@uga.edu)
Plant Genome Mapping Laboratory, University of Georgia

This script:
  1) Checks for the hymenoptera_odb10 (or user-chosen) BUSCO lineage dataset. 
     If not found, attempts to download it.
  2) Reads a TSV of genomes (genome_list.tsv).
  3) For each genome:
     - Runs BUSCO (genome mode).
     - Runs DIAMOND alignment of predicted proteins vs. D. melanogaster known proteins.
     - Calls Count_pident.py to summarize pident distributions.
     - Collects BUSCO completeness into a master_summary.tsv.

Parallel usage:
  - Set --cpu for multi-threading on each genome 
    (BUSCO and DIAMOND will use that many threads).
  - For HPC clusters with SLURM, see the example sbatch script below.

Usage (local run):
  python beeome_minimal_pipeline.py --genome_list genome_list.tsv \
                                    --lineage hymenoptera_odb10 \
                                    --dmel_faa dmel_known_refseq.faa \
                                    --cpu 8
"""

import os
import sys
import subprocess
import argparse
import logging
import re
from datetime import datetime


def setup_logging(logfile=None):
    """
    Configures logging to both console and an optional logfile.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", 
                                  datefmt="%Y-%m-%d %H:%M:%S")

    # Console output
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if logfile:
        # File output
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def check_and_download_lineage(lineage_name):
    """
    Checks if the lineage folder (e.g., 'hymenoptera_odb10') exists locally.
    If not, tries to download it with 'busco download --lineages <lineage_name>'.
    Adjust as needed for BUSCO versions.
    """
    if os.path.isdir(lineage_name):
        logging.info(f"BUSCO lineage folder '{lineage_name}' found. No download needed.")
        return lineage_name

    # If not found, attempt download
    logging.info(f"BUSCO lineage folder '{lineage_name}' not found. Attempting to download...")
    try:
        cmd = ["busco", "download", "--lineages", lineage_name]
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        logging.error("Failed to download the BUSCO lineage data. Check internet or lineage spelling.")
        sys.exit(1)

    if not os.path.isdir(lineage_name):
        logging.error(f"Download completed, but the folder '{lineage_name}' is still not found.")
        sys.exit(1)

    return lineage_name


def run_busco(genome_fasta, output_prefix, lineage_path, cpu=4):
    """
    Runs BUSCO in 'genome' mode using the specified lineage_path.
    """
    cmd = [
        "busco",
        f"--in={genome_fasta}",
        f"--out={output_prefix}",
        f"--lineage_path={lineage_path}",
        "--mode=genome",
        f"--cpu={cpu}"
    ]
    logging.info("Running BUSCO with command:\n  " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_busco_summary(prefix):
    """
    Parses the short_summary_<prefix>.txt from run_<prefix>/ to extract stats.
    Returns a dict of:
      Complete(%), Single(%), Duplicated(%), Fragmented(%), Missing(%), TotalBUSCO
    """
    summary_file = f"run_{prefix}/short_summary_{prefix}.txt"
    stats = {
        "Complete(%)": 0.0,
        "Single(%)": 0.0,
        "Duplicated(%)": 0.0,
        "Fragmented(%)": 0.0,
        "Missing(%)": 0.0,
        "TotalBUSCO": 0
    }

    if not os.path.exists(summary_file):
        logging.warning(f"BUSCO summary file not found: {summary_file}")
        return stats

    pattern = r"C:([\d\.]+)%\[S:([\d\.]+)%,D:([\d\.]+)%\],F:([\d\.]+)%,M:([\d\.]+)%,n:(\d+)"
    with open(summary_file, "r") as f:
        for line in f:
            if "C:" in line and "S:" in line and "D:" in line and "n:" in line:
                match = re.search(pattern, line)
                if match:
                    stats["Complete(%)"]   = float(match.group(1))
                    stats["Single(%)"]     = float(match.group(2))
                    stats["Duplicated(%)"] = float(match.group(3))
                    stats["Fragmented(%)"] = float(match.group(4))
                    stats["Missing(%)"]    = float(match.group(5))
                    stats["TotalBUSCO"]    = int(match.group(6))
    return stats


def diamond_align(dmel_faa, query_faa, out_tsv, cpu=4):
    """
    Uses DIAMOND to align query_faa (predicted proteins) 
    against dmel_faa (D. melanogaster known proteins).
    """
    # Step 1: Make the DIAMOND database 
    cmd_db = [
        "diamond", "makedb",
        "--in", dmel_faa,
        "-d", "dmel_db"
    ]
    logging.info("Creating DIAMOND DB:\n  " + " ".join(cmd_db))
    subprocess.run(cmd_db, check=True)

    # Step 2: DIAMOND blastp
    cmd_blast = [
        "diamond", "blastp",
        "--db", "dmel_db.dmnd",
        "--query", query_faa,
        "--out", out_tsv,
        "--outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "--threads", str(cpu)
    ]
    logging.info("DIAMOND alignment:\n  " + " ".join(cmd_blast))
    subprocess.run(cmd_blast, check=True)


def run_count_pident(blast_tsv, out_prefix):
    """
    Calls Count_pident.py on the BLAST/DIAMOND TSV to generate 
    histogram and cumulative distribution of pident.
    """
    script = "./Count_pident.py"
    if not os.path.exists(script):
        logging.warning(f"Count_pident.py not found at {script}. Skipping pident analysis.")
        return
    cmd = ["python", script, blast_tsv, out_prefix]
    logging.info("Running Count_pident.py:\n  " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Minimal Beeome pipeline, with parallel usage & BUSCO lineage check. By Ankush Sharma (UGA)."
    )
    parser.add_argument("--genome_list", required=True,
                        help="TSV with: GenomeID GenomeFasta ProteinsFaa [OutputPrefix].")
    parser.add_argument("--lineage", default="hymenoptera_odb10",
                        help="Name (or path) of the BUSCO lineage data to use.")
    parser.add_argument("--dmel_faa", default="dmel_known_refseq.faa",
                        help="Path to D. melanogaster known proteins FASTA.")
    parser.add_argument("--cpu", type=int, default=4,
                        help="Number of threads for BUSCO/DIAMOND.")
    parser.add_argument("--log", default="minimal_pipeline.log",
                        help="Log file name.")
    args = parser.parse_args()

    logger = setup_logging(args.log)
    start_time = datetime.now()

    # 1) Check or download the lineage dataset
    lineage_path = check_and_download_lineage(args.lineage)

    # 2) Prepare master summary
    out_summary = open("master_summary.tsv", "w")
    out_summary.write("GenomeID\tComplete(%)\tSingle(%)\tDuplicated(%)\tFragmented(%)\tMissing(%)\tTotalBUSCO\n")

    # 3) Read the genome list
    if not os.path.exists(args.genome_list):
        logger.error(f"Genome list file not found: {args.genome_list}")
        sys.exit(1)

    with open(args.genome_list, "r") as fh:
        for line in fh:
            if not line.strip() or line.strip().startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                logger.warning("Skipping line with too few columns:\n" + line)
                continue

            genome_id   = parts[0]
            genome_fasta= parts[1]
            proteins_faa= parts[2]
            if len(parts) > 3:
                busco_prefix = parts[3]
            else:
                busco_prefix = genome_id

            logger.info(f"\nProcessing genome: {genome_id}")
            logger.info(f"  Genome FASTA = {genome_fasta}")
            logger.info(f"  Proteins FASTA = {proteins_faa}")
            logger.info(f"  BUSCO prefix = {busco_prefix}")

            # --- BUSCO ---
            run_busco(genome_fasta, busco_prefix, lineage_path=lineage_path, cpu=args.cpu)

            # --- Parse BUSCO results ---
            stats = parse_busco_summary(busco_prefix)

            # --- DIAMOND alignment ---
            if os.path.exists(proteins_faa):
                diamond_out = f"{busco_prefix}_diamond.tsv"
                diamond_align(args.dmel_faa, proteins_faa, diamond_out, cpu=args.cpu)
                # --- pident distribution ---
                run_count_pident(diamond_out, f"{busco_prefix}_pident")
            else:
                logger.warning(f"Proteins file not found: {proteins_faa}")

            # Write stats
            out_summary.write(
                f"{genome_id}\t{stats['Complete(%)']}\t{stats['Single(%)']}\t{stats['Duplicated(%)']}\t"
                f"{stats['Fragmented(%)']}\t{stats['Missing(%)']}\t{stats['TotalBUSCO']}\n"
            )

    out_summary.close()
    elapsed = datetime.now() - start_time
    logger.info(f"\nAll genomes processed! Total pipeline time: {elapsed}\n")


if __name__ == "__main__":
    main()
