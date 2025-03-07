
# Genome QC and Contamination Checking

This repository contains **two major pipelines** for genome quality control:

1. **`busco_diamond.py`** (a single Python script)  
   - Runs [BUSCO](https://busco.ezlab.org/) to assess completeness of each genome.  
   - (Optionally) runs **DIAMOND** alignment of predicted proteins against *D. melanogaster* known RefSeq proteins and calls `Count_pident.py` for percent-identity analysis.

2. **Contamination Pipeline** (`diamond_contamination.py` + `run_contamination.sbatch`)  
   - Screens genomes for known contaminants (e.g., UniVec or PhiX proteins) using **DIAMOND** in *blastx* style.  
   - Merges results into a single file, labeling each genome and contaminant source.



