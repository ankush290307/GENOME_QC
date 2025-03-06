#!/bin/bash
#SBATCH --job-name=beeome_pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=beeome_pipeline_%j.log
# Add parameters basd on the HPC used
# Load modules (Module spider python busco diamond) then load the specific package below
# module load python
# module load busco
# module load diamond


echo "Starting pipeline job on $(hostname) at $(date)"

python beeome_minimal_pipeline.py \
    --genome_list genome_list.tsv \
    --lineage hymenoptera_odb10 \
    --dmel_faa dmel_known_refseq.faa \
    --cpu $SLURM_CPUS_PER_TASK \
    --log "beeome_pipeline.log"

echo "Pipeline completed at $(date)"
