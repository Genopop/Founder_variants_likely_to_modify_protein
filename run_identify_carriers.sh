#!/bin/bash
#SBATCH --account=<your_account>
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=5G
#SBATCH --time=0-01:00

for chr in {1..22}; do
echo "=== Chromosome $chr ==="
julia identify_carriers.jl /path/to/wgs_cluster1_rfd10_5carriers_CCDS_flanking_maf10_chr${chr}.tped $OUT_DIR/wgs_cluster1_rfd10_5carriers_CCDS_flanking_maf10_chr${chr}.tfam $OUT_DIR/snps_rfd10_chr${chr}.txt $OUT_DIR/carriers_by_variant_snps_rfd10_5carriers_CCDS_flanking_maf10_chr${chr}.txt
done
