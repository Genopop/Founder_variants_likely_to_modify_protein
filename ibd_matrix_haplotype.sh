#!/bin/bash
#SBATCH --account=<your_account>
#SBATCH --mem=100G
#SBATCH --time=0-01:00
#SBATCH --cpus-per-task=10
#SBATCH --array=1-22
#SBATCH --job-name=ibd_chr

CHR=${SLURM_ARRAY_TASK_ID}

BIM_FILE="/path/to/.bim"
TPED_FILE="/path/to/chr${CHR}.tped"
TFAM_FILE="/path/to/chr${CHR}.tfam"
IBD_FILE="/path/to/refinedibd_chr${CHR}.ibd.gz"

julia ibd_select_one_haplotype_per_carrier.jl "$CHR" "$TPED_FILE" "$TFAM_FILE" "$IBD_FILE"
