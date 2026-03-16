### Here are usage examples. For detailed parameters, please refer to the source code.

#01 filter
grep -v "inter" new_ins.bed | awk '{print $5"\t"$6"\t"$7}' > final_hosa_ins.bed
grep -v "inter" new_del.bed | awk '{print $1"\t"$2"\t"$3}' > final_hosa_del.bed

##02 Extract the coordinate blocks (only the forward strand) of the reference species and the target species for comparison
python extract_maf_bed.py -i filter_realign.maf -s mu_mu -o ref_mu_mu_block.bed -f bed
python extract_maf_bed.py -i filter_realign.maf -s ho_sa -o target_ho_sa_block.bed -f bed

##03 Extract the final SV of the detected SV in the alignment block
intersectBed -a final_hosa_ins.bed -b target_ho_sa_block.bed -wa -f 1 -u > final_human_ins.bed
intersectBed -a final_hosa_del.bed -b ref_mu_mu_block.bed -wa -f 1 -u > final_human_del.bed

##04 Obtain the corresponding variant sequences of the reference species and the target species
python obtain_SV_sequence.py --bed final_human_del.bed --maf filter_realign.maf --output del --ref-species mu_mu --target-species ho_sa --extend 10 --threshold 0.5
python obtain_SV_sequence.py --bed final_human_ins.bed --maf filter_realign.maf --output ins --ref-species ho_sa --target-species mu_mu --extend 10 --threshold 0.5

##05 Extract the block where SV is located
python check_SV_obtain_MAFblock.py --bed final_human_del.bed --maf filter_realign.maf --output block_check_del.txt --context 1 --species mu_mu --search-range 1000
python check_SV_obtain_MAFblock.py --bed final_human_ins.bed --maf filter_realign.maf --output block_check_ins.txt --context 1 --species ho_sa --search-range 1000