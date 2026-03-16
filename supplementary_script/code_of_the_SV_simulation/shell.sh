### Here are usage examples. For detailed parameters, please refer to the source code. Among them, vital_SV.bed contains the sequence editing position information that needs to be prepared in advance; vital_100_SV_simulate.bed contains the sequence insertion position information that needs to be prepared in advance. Additionally, it is necessary to pay attention to the path of the output information.

## 01 Generate 2500 homologous sequences with a sequence difference of 0 for 2k and 20k.
cd ./identity_100_seq
python generate_holomogy_sequences.py -l 1000 -d 0.0 -n 1250 -o homologous_sequences_1k_100.fa
python generate_holomogy_sequences.py -l 20000 -d 0.0 -n 1250 -o homologous_sequences_2k_100.fa

## 02 Edit the sequences in the homologous sequences, and prepare sequence files with INS, DEL, and INV (2 sequences); sequence files with TRA (3 sequences)
cd ..
python /work/home/acm2g39zj0/wangyj/SV/new_SV_method/00_script/generate_sv.py -b ./vital_SV.bed -i ./identity_100_seq -o ./identity_100_seq_SV

## 03 Generate a simulated genome. The target genome needs to be modified separately. Open the "generate_genome.py" file and increase the number of chromosomes by one.
python generate_genome.py -o genome1.fa -s 12850000
python generate_genome.py -o genome2.fa -s 12850000
python generate_genome.py -o target.fa -s 12850000
python generate_genome.py -o genome3.fa -s 12850000

## 04 Insert the homologous sequence file containing SV into the simulated genome
python insert_seq_to_genome.py --aln-dir ./identity_100_seq_SV --target-dir ./identity_100_genome --bed-file ./vital_100_SV_simulate.bed --output-dir ./final_genome --coord-file vital_100_SV_final.bed

## 05 Extract the reference species coordinates of the simulated SV
cd ..
python extract_SV_simulation_bed.py vital_100_SV_marker.bed ./vital_100_SV_final.bed vital_100_SV_final_check.bed
grep -i "genome1" vital_100_SV_final_check.bed | grep -i "DEL" > vital_100_del.bed
grep -i "genome1" vital_100_SV_final_check.bed | grep -i "INS" > vital_100_ins.bed;cat vital_100_ins.bed | awk '{print $1"\t"$2"\t"$3+$5"\t"$4"\t"$5}' > vital_100_ins_gai.bed
grep -i "genome1" vital_100_SV_final_check.bed | grep -i "intraTRA" > vital_100_intraTRA.bed;cat vital_100_intraTRA.bed | awk '{print $1"\t"$2"\t"$3+$5"\t"$4"\t"$5}' > vital_100_intraTRA_vary.bed
grep -i "genome1" vital_100_SV_final_check.bed | grep -i "interTRA" > vital_100_interTRA.bed;cat vital_100_interTRA.bed | awk '{print $1"\t"$2"\t"$3+$5"\t"$4"\t"$5}' > vital_100_interTRA_vary.bed
grep -i "genome1" vital_100_SV_final_check.bed | grep -i "INV" > vital_100_inv.bed