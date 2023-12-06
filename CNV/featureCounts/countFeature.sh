#!/bin/bash
#SBATCH --job-name=countGene
#SBATCH --time=1:00:00
#SBATCH --output=countGene_%A_%a.out
#SBATCH --error=countGene_%A_%a.out
#SBATCH --array=2-282
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000
#SBATCH --account=heqixin
#SBATCH --mail-type=ALL
#SBATCH --mail-user=heqixin@purdue.edu

# Set the line number you want to read
line_number=$SLURM_ARRAY_TASK_ID  
# Specify the file path
file_path="bam_list.txt"
# Use sed to extract the desired line
desired_sample=$(sed -n "${line_number}p" "$file_path")

./bin/featureCounts -a PlasmoDB-65_Pfalciparum3D7.gff -o counts_${SLURM_ARRAY_TASK_ID} -p -B -C $desired_sample
