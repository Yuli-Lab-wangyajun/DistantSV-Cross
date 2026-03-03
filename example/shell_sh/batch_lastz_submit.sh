#!/bin/bash
# Batch submission comparison job script

echo "Start submitting tasks..."

echo "submitting tasks: genome1_chr1_genome2_lastz.sh"
sbatch genome1_chr1_genome2_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr1_genome3_lastz.sh"
sbatch genome1_chr1_genome3_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr1_target_lastz.sh"
sbatch genome1_chr1_target_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr2_genome2_lastz.sh"
sbatch genome1_chr2_genome2_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr2_genome3_lastz.sh"
sbatch genome1_chr2_genome3_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr2_target_lastz.sh"
sbatch genome1_chr2_target_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr3_genome2_lastz.sh"
sbatch genome1_chr3_genome2_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr3_genome3_lastz.sh"
sbatch genome1_chr3_genome3_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr3_target_lastz.sh"
sbatch genome1_chr3_target_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr4_genome2_lastz.sh"
sbatch genome1_chr4_genome2_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr4_genome3_lastz.sh"
sbatch genome1_chr4_genome3_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr4_target_lastz.sh"
sbatch genome1_chr4_target_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr5_genome2_lastz.sh"
sbatch genome1_chr5_genome2_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr5_genome3_lastz.sh"
sbatch genome1_chr5_genome3_lastz.sh
sleep 1
echo "submitting tasks: genome1_chr5_target_lastz.sh"
sbatch genome1_chr5_target_lastz.sh
sleep 1

