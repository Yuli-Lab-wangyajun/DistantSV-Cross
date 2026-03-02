#!/bin/bash
#SBATCH -J JOB_NAME_PLACEHOLDER
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o /work/home/acm2g39zj0/err_out/std.out.%j
#SBATCH -e /work/home/acm2g39zj0/err_out/std.err.%j

source ~/miniconda3/bin/activate
conda activate test
cd WORK_DIR_PLACEHOLDER/axt/CHROMOSOME_PLACEHOLDER_axt
lastz WORK_DIR_PLACEHOLDER/ref_chr/CHROMOSOME_PLACEHOLDER GENOME_NAME_PLACEHOLDER --notransition --ydrop=15000 --gappedthresh=3000 --hspthresh=4500 --format=axt --out=GENOME_NAME_OUTPUT.axt