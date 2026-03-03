#!/bin/bash
# Batch submission script for SV analysis jobs

echo "Starting to submit SV analysis jobs..."

echo "Submitting job: genome1_chr1_SV.sh"
sbatch genome1_chr1_SV.sh
sleep 2
echo "Submitting job: genome1_chr2_SV.sh"
sbatch genome1_chr2_SV.sh
sleep 2
echo "Submitting job: genome1_chr3_SV.sh"
sbatch genome1_chr3_SV.sh
sleep 2
echo "Submitting job: genome1_chr4_SV.sh"
sbatch genome1_chr4_SV.sh
sleep 2
echo "Submitting job: genome1_chr5_SV.sh"
sbatch genome1_chr5_SV.sh
sleep 2

echo "All SV analysis jobs have been submitted!"
