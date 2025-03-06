#!/bin/bash
#SBATCH --job-name=diamond_contam
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=diamond_contamination_%j.log

################################################################################
# Author: Ankush Sharma (Ankush.Sharma@uga.edu)
# This SLURM script runs "diamond_contamination.py" twice: 
#   1) For UniVec proteins
#   2) For PhiX proteins
# and merges both result files into a final table.
################################################################################

# 1) Load or activate DIAMOND use module spider diamond or conda env
module load diamond

# 2) Variables for file paths
GENOME_LIST="genome_list.tsv"         # The two-column file: GenomeID <tab> /path/to/genome.fna
UNIVEC_FASTA="univec_proteins.faa"    # Hypothetical protein version of UniVec
PHIX_FASTA="phix_proteins.faa"        # Protein reference for PhiX
THREADS=$SLURM_CPUS_PER_TASK

# 3) Run diamond_contamination.py for UniVec
python diamond_contamination.py \
  --genome_list "$GENOME_LIST" \
  --ref_faa "$UNIVEC_FASTA" \
  --contaminant_label "univec" \
  --out univec_compiled.tsv \
  --threads "$THREADS"

# 4) Run diamond_contamination.py for PhiX
python diamond_contamination.py \
  --genome_list "$GENOME_LIST" \
  --ref_faa "$PHIX_FASTA" \
  --contaminant_label "phix" \
  --out phix_compiled.tsv \
  --threads "$THREADS"

# 5) Combine both results into one final file:
#    We'll copy the header line from the first file, then skip the header from the second.
#    Or we can do an AWK approach that prints the header from the first file only once.

FINAL_FILE="final_combined.tsv"

# Grab the header from univec_compiled.tsv
head -n 1 univec_compiled.tsv > "$FINAL_FILE"

# Append data from univec_compiled.tsv (excluding header line) 
tail -n +2 univec_compiled.tsv >> "$FINAL_FILE"

# Append data from phix_compiled.tsv (excluding header line)
tail -n +2 phix_compiled.tsv >> "$FINAL_FILE"

echo "All done! Final combined file is $FINAL_FILE"

################################################################################
# End of script
################################################################################
