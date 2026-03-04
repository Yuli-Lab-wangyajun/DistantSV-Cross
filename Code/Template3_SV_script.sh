#!/bin/bash
#SBATCH -J JOB_NAME_PLACEHOLDER
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o your_process_file_directory/JOB_NAME_PLACEHOLDER.out
#SBATCH -e your_process_file_directory/JOB_NAME_PLACEHOLDER.err

conda activate your_env

cd WORK_DIR_PLACEHOLDER/CHROMOSOME_PLACEHOLDER_multiz_out/multiz_out
python ../../Module1_MAF_filter.py --input CHROMOSOME_PLACEHOLDER.maf --output CHROMOSOME_PLACEHOLDER_filter.maf --min_sequences 4 --ref_species REF_GENOME_PLACEHOLDER --min_ref_length 50 --max_gap_ratio 0.1 --gap_region_length 100
python ../../Module2_MAF_realigner.py -i CHROMOSOME_PLACEHOLDER_filter.maf -o CHROMOSOME_PLACEHOLDER_filter_realign.maf
python ../../Module3_call_SV.py --maf CHROMOSOME_PLACEHOLDER_filter_realign.maf --target TARGET_SPECIES_PLACEHOLDER --reference REF_GENOME_PLACEHOLDER --non_target "NON_TARGET_PLACEHOLDER" --output_prefix CHROMOSOME_PLACEHOLDER --max_target_gap 10 --min_ref_gap 50 --min_target_gap 50 --max_ref_gap 10 --max_gap 100 --inv_more_block_max 1000 --min_intra_size 50 --max_gap_ratio 0.3
python ../../Module4_sequence_consistency_calculation.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_insertions.bed -o CHROMOSOME_PLACEHOLDER_ins.bed -r REF_GENOME_PLACEHOLDER
python ../../Module4_sequence_consistency_calculation.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_deletions.bed -o CHROMOSOME_PLACEHOLDER_del.bed -r REF_GENOME_PLACEHOLDER
python ../../Module4_sequence_consistency_calculation.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_inversions.bed -o CHROMOSOME_PLACEHOLDER_inv.bed -r REF_GENOME_PLACEHOLDER
python ../../Module4_sequence_consistency_calculation.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_interchrom_translocations.bed -o CHROMOSOME_PLACEHOLDER_interTRA.bed -r REF_GENOME_PLACEHOLDER
python ../../Module4_sequence_consistency_calculation.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_intrachrom_translocations.bed -o CHROMOSOME_PLACEHOLDER_intraTRA.bed -r REF_GENOME_PLACEHOLDER

