#!/bin/bash 
#SBATCH -J JOB_NAME_PLACEHOLDER
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o /work/home/acm2g39zj0/err_out/multiz_CHROMOSOME_PLACEHOLDER.out
#SBATCH -e /work/home/acm2g39zj0/err_out/multiz_CHROMOSOME_PLACEHOLDER.err

source ~/miniconda3/bin/activate
conda activate test
cd WORK_DIR_PLACEHOLDER
nib=WORK_DIR_PLACEHOLDER/species_nibs
sizes=WORK_DIR_PLACEHOLDER/species_sizes
bit=WORK_DIR_PLACEHOLDER/species_2bit
chr_axt=WORK_DIR_PLACEHOLDER/axt/CHROMOSOME_PLACEHOLDER_axt
axt_list=`ls $chr_axt`

mkdir -p CHROMOSOME_PLACEHOLDER_multiz_out
cd CHROMOSOME_PLACEHOLDER_multiz_out

for i in $axt_list
do
axtChain -linearGap=loose -minScore=3000 $chr_axt/$i $nib/REF_GENOME_PLACEHOLDER.nib $nib/${i%.axt}.nib ${i%.axt}.chain
chainPreNet ${i%.axt}.chain $sizes/REF_GENOME_PLACEHOLDER.sizes $sizes/${i%.axt}.sizes ${i%.axt}.pre.chain
chainNet ${i%.axt}.pre.chain $sizes/REF_GENOME_PLACEHOLDER.sizes $sizes/${i%.axt}.sizes stdout /dev/null | netSyntenic stdin ${i%.axt}.chain.net
netToAxt ${i%.axt}.chain.net ${i%.axt}.pre.chain $bit/REF_GENOME_PLACEHOLDER.2bit $bit/${i%.axt}.2bit ${i%.axt}.Axt
axtToMaf ${i%.axt}.Axt $sizes/REF_GENOME_PLACEHOLDER.sizes $sizes/${i%.axt}.sizes REF_GENOME_PLACEHOLDER.${i%.axt}.sing.maf -tPrefix=REF_GENOME_PLACEHOLDER. -qPrefix=${i%.axt}.
done

mkdir -p multiz_out
roast + E=REF_GENOME_PLACEHOLDER T=WORK_DIR_PLACEHOLDER/CHROMOSOME_PLACEHOLDER_multiz_out/multiz_out 'SPECIES_TREE_PLACEHOLDER' WORK_DIR_PLACEHOLDER/CHROMOSOME_PLACEHOLDER_multiz_out/*.sing.maf WORK_DIR_PLACEHOLDER/CHROMOSOME_PLACEHOLDER_multiz_out/multiz_out/CHROMOSOME_PLACEHOLDER.maf