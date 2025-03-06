# BUSCO & Contamination Screening Pipelines

This repository contains two pipelines for genome quality control:

1. **BUSCO Pipeline** – Runs [BUSCO](https://busco.ezlab.org/) on a list of genomes to assess completeness and compiles stats in a master summary. It can also (optionally) run **DIAMOND** on predicted proteins against _Drosophila melanogaster_ known RefSeq proteins to examine percent identity distribution.
2. **DIAMOND Contamination Pipeline** – Screens genome assemblies for contaminant proteins (e.g., UniVec or PhiX) by translating the genome (blastx-style) and comparing to a contaminant protein reference.

---

## Files

- **`beeome_minimal_pipeline.py`**  
  - Main BUSCO pipeline script.
  - Reads a genome list, runs BUSCO, and (if provided) runs **DIAMOND** alignment against _Drosophila melanogaster_ known proteins for identity checks.
  - Produces `master_summary.tsv` (BUSCO stats) + optional DIAMOND alignment outputs.

- **`diamond_contamination.py`**  
  - Screens genomes for contamination by running DIAMOND blastx against a *protein* reference (e.g., UniVec or PhiX).
  - Merges all hits into a single output, prefixing lines with the `GenomeID` (and optionally a label).

- **`run_contamination.sbatch`**  
  - Example SLURM script that calls `diamond_contamination.py` **twice**, once for UniVec (or some vector reference) and once for PhiX. Then merges both results.

- **`genome_list.tsv`**  
  - Example two-column TSV with:
    ```
    #GenomeID    /path/to/genome.fna
    sampleA      /scratch/user/sampleA.fna
    sampleB      /scratch/user/sampleB.fna
    ```
  - Lines beginning with `#` are ignored.

---

## Requirements

- **Python 3**  
- **BUSCO** (if using the BUSCO pipeline)  
- **DIAMOND** (if performing alignment steps against _Drosophila melanogaster_ or contaminant references)  
- **SLURM** (for HPC usage)

---

## Usage

### 1) BUSCO Pipeline (with optional DIAMOND on D. melanogaster)

1. **Genome List**: create `genome_list.tsv` with lines like:
