#!/bin/bash
#SBATCH -J JOB_NAME_PLACEHOLDER
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o your_process_file_directory/JOB_NAME_PLACEHOLDER.out
#SBATCH -e your_process_file_directory/JOB_NAME_PLACEHOLDER.err

conda activate your_env

cd WORK_DIR_PLACEHOLDER/axt/CHROMOSOME_PLACEHOLDER_axt
lastz WORK_DIR_PLACEHOLDER/ref_chr/CHROMOSOME_PLACEHOLDER GENOME_NAME_PLACEHOLDER --notransition --ydrop=15000 --gappedthresh=3000 --hspthresh=4500 --format=axt --out=GENOME_NAME_OUTPUT.axt

