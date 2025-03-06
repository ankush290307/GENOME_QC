#!/bin/bash
#SBATCH --job-name=beeome_QC_pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=beeome_QC_pipeline_%j.log
# Add the HPC specifcic sbatch parameters above
# Load the required modules (Do module spider python busco diamond ) and then load the latest version available below
# module load python
# module load busco
# module load diamond



echo "Starting pipeline job on $(hostname) at $(date)"

python busco_diamond.py \
    --genome_list genome_list.tsv \
    --lineage hymenoptera_odb10 \
    --dmel_faa dmel_known_refseq.faa \
    --cpu $SLURM_CPUS_PER_TASK \
    --log "beeome_pipeline.log"

echo "Pipeline completed at $(date)"
