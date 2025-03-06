
# BUSCO & Contamination Screening Pipelines

This repository contains two pipelines for genome quality control:

1. **BUSCO Pipeline** – Runs BUSCO on a list of genomes and compiles completeness stats.
2. **DIAMOND Contamination Pipeline** – Screens genome assemblies for contaminants (e.g., UniVec and PhiX) using DIAMOND blastx and merges results.

---

## Files

- **`beeome_minimal_pipeline.py`**  
  Runs BUSCO on genomes (requires a TSV with GenomeID and genome FASTA paths).

- **`diamond_contamination.py`**  
  Screens genomes for contamination by running DIAMOND blastx against a contaminant protein reference.  
  Accepts an extra `--contaminant_label` parameter to tag results.

- **`run_contamination.sbatch`**  
  SLURM submission script that runs the contamination screening twice (once for UniVec and once for PhiX) and then merges the outputs.

- **`genome_list.tsv`**  
  Example TSV file with two columns (tab-delimited):
