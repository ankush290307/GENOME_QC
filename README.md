
# Genome QC and Contamination Checking

This repository contains **two major pipelines** for genome quality control:

1. **`busco_diamond.py`** (a single Python script)  
   - Runs [BUSCO](https://busco.ezlab.org/) to assess completeness of each genome.  
   - (Optionally) runs **DIAMOND** alignment of predicted proteins against *D. melanogaster* known RefSeq proteins and calls `Count_pident.py` for percent-identity analysis.

2. **Contamination Pipeline** (`diamond_contamination.py` + `run_contamination.sbatch`)  
   - Screens genomes for known contaminants (e.g., UniVec or PhiX proteins) using **DIAMOND** in *blastx* style.  
   - Merges results into a single file, labeling each genome and contaminant source.

Below are details on each pipeline.

---

## 1) `busco_diamond.py`

### Overview

- **Checks** for the BUSCO lineage dataset (default `hymenoptera_odb10`), downloads if missing.
- **Reads** a TSV (`genome_list.tsv`) listing:
#GenomeID /path/to/genome.fna /path/to/predicted_proteins.faa [OutputPrefix?]
- **For each genome**:
1. Runs **BUSCO** in genome mode, storing results in `run_<prefix>/`.
2. If predicted proteins are provided, runs **DIAMOND** alignment vs. a *Drosophila melanogaster* reference (`--dmel_faa`).
3. Calls `Count_pident.py` (if present) to produce histograms/cumulative distributions of alignment identity.
4. Collects BUSCO completeness stats into `master_summary.tsv`.

### Requirements
- **Python 3**
- **BUSCO** (and a downloaded lineage, e.g. `hymenoptera_odb10`)
- **DIAMOND** (if doing the alignment steps)
- (Optional) `Count_pident.py` to parse alignment results


