#!/bin/bash
#SBATCH -J script_sv
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o /work/home/acm2g39zj0/err_out/std.out.%j
#SBATCH -e /work/home/acm2g39zj0/err_out/std.err.%j

source ~/miniconda3/bin/activate
conda activate test
##step1 Genome alignment preprocessing
cd /work/home/acm2g39zj0/wangyj/SV/new_SV_method/really_data/script_test
dos2unix *.sh *.py
./01_pre.sh -s species.txt -r genome1 -d /work/home/acm2g39zj0/wangyj/SV/new_SV_method/really_data/script_test

##step2 Pairwise genome alignment
./02_align.sh -s species.txt -r genome1 -d /work/home/acm2g39zj0/wangyj/SV/new_SV_method/really_data/script_test -t Template1_lastz_script.sh

##step3 Generating progressive Multi-Sequence Alignment Files
./03_MSA.sh -s species.txt -r genome1 -d /work/home/acm2g39zj0/wangyj/SV/new_SV_method/really_data/script_test -t Template2_multiz_script.sh -tree "(((genome1 genome2) target) genome3)"

##step4-6 Alignment block filtering and realignment; Calling structural variation; Calculation of the Mean Sequence Consistency
./04_SV.sh -s species.txt -r genome1 -d /work/home/acm2g39zj0/wangyj/SV/new_SV_method/really_data/script_test -t Template3_SV_script.sh -target target -non_target "genome1,genome2,genome3"

