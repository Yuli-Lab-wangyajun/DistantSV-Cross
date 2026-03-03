#!/bin/bash
# Batch submission script for multiple sequence alignment jobs

echo "Starting to submit multiple sequence alignment jobs..."

echo "Submitting job: genome1_chr1_multiz.sh"
sbatch genome1_chr1_multiz.sh
sleep 2
echo "Submitting job: genome1_chr2_multiz.sh"
sbatch genome1_chr2_multiz.sh
sleep 2
echo "Submitting job: genome1_chr3_multiz.sh"
sbatch genome1_chr3_multiz.sh
sleep 2
echo "Submitting job: genome1_chr4_multiz.sh"
sbatch genome1_chr4_multiz.sh
sleep 2
echo "Submitting job: genome1_chr5_multiz.sh"
sbatch genome1_chr5_multiz.sh
sleep 2

echo "All multiple sequence alignment jobs have been submitted!"
