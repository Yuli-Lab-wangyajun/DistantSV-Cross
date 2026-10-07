#!/bin/bash
#SBATCH -J control_script.sh
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o your_process_file_directory/control_script.out
#SBATCH -e your_process_file_directory/control_script.err

conda activate your_env

cd working_directory
dos2unix *.sh *.py
##step1 Genome alignment preprocessing
./01_pre.sh -s species.txt -r ref -d working_directory

##step2 Pairwise genome alignment
./02_align.sh -s species.txt -r ref -d working_directory -t Template1_lastz_script.sh

##step3 Generating progressive Multi-Sequence Alignment Files
./03_MSA.sh -s species.txt -r ref -d working_directory -t Template2_multiz_script.sh -tree "(((ref query1) target) query2)"

##step4-6 Alignment block filtering and realignment; Calling structural variation; Calculation of the Mean Sequence Consistency
./04_SV.sh -s species.txt -r ref -d working_directory -t Template3_SV_script.sh -target target -non_target "ref,query1,query2"

python GNUparallel.py working_directory/shell_sh -s ref,query1,target,query2 -r ref -g 60