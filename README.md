## Overview
DistantSV-Cross is a computational pipeline specifically designed for detecting structural variations among distantly related multiple species. It primarily comprises six computational steps to identify species-specific genomic structural variations in a target species relative to other species, ultimately generating insertion, deletion, inversion, intra-chromosomal translocation, and inter-chromosomal translocation information in BED format. Compared to other variant detection tools based on pairwise species comparison via genome assembly, as well as existing multi-species structural variation detection software, DistantSV-Cross demonstrates superior performance, achieving higher precision, recall, and F1-score. <br/>

Table of Contents:
- [Overview](#overview)
- [Installation:](#installation)
  - [1. Download the source code](#1-download-the-source-code)
  - [2. Install the required dependencies for running the software](#2-install-the-required-dependencies-for-running-the-software)
- [Data preparation](#data-preparation)
  - [1. Prepare the species.txt file](#1-prepare-the-speciestxt-file)
  - [2. Determine the phylogenetic relationships among the species](#2-determine-the-phylogenetic-relationships-among-the-species)
- [Usage](#usage)
  - [1. Genome alignment preprocessing](#1-genome-alignment-preprocessing)
  - [2. Pairwise genome alignment](#2-pairwise-genome-alignment)
  - [3. Generating Progressive Multi-Sequence Alignment Files](#3-generating-progressive-multi-sequence-alignment-files)
  - [4-5-6. Alignment block filtering and realignment](#4-5-6-alignment-block-filtering-and-realignment)
- [Detailed parameters](#detailed-parameters)
- [License](#license)

## Installation:
DistantSV-Cross provides a complete pipeline for running the code, requiring only the installation of the required software, Python packages, and downloading the source code provided by us.

### 1. Download the source code
```
wget https://github.com/Yuli-Lab-wangyajun/DistantSV-Cross.git
or
git clone https://github.com/Yuli-Lab-wangyajun/DistantSV-Cross
```

### 2. Install the required dependencies for running the software
```
## It is recommended to create a new environment to install the required software or packages. The Python environment should be 3.8+, here we use env_test as an example.
conda create -n env_test python=3.8 pandas numpy -y

## Manually install the following software and add them to the env_test environment:
lastz, multiz, muscle, mafft

## From the installation packages at (https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/), find the following packages and add them to the env_test environment:
axtChain, chainPreNet, chainNet, netToAxt, axtToMaf, faSize, faToTwoBit, faSplit
```

## Data preparation
I will use the sample files in the demo folder to demonstrate the usage steps of DistantSV-Cross in detail. All the following steps should be performed in the working directory "/Yuli-Lab-wangyajun/DistantSV-Cross/example".

Before running, the following information needs to be prepared:
### 1. Prepare the species.txt file
```
genome1	/Yuli-Lab-wangyajun/DistantSV-Cross/example/genomes_path/genome1.fa
genome2	/Yuli-Lab-wangyajun/DistantSV-Cross/example/genomes_path/genome2.fa
genome3	/Yuli-Lab-wangyajun/DistantSV-Cross/example/genomes_path/genome3.fa
target	/Yuli-Lab-wangyajun/DistantSV-Cross/example/genomes_path/target.fa
```
Note: The first column is the abbreviated name for each species' genome, which you can customize for your convenience—it is recommended to keep it concise. The second column is the absolute path to the genome. Please note that these two columns must be separated by a tab. In this example, target.fa is designated as the target species to be analyzed, and the reference species is genome1.fa.

### 2. Determine the phylogenetic relationships among the species
```
(((genome1 target) genome2) genome3);
```
Note: The phylogenetic relationships in Newick format are primarily used for constructing multiple sequence alignments in subsequent steps with multiz.

## Usage 
After preparing the above information, SV detection can be formally performed. In this example, the reference genome is genome1.fa, and the target genome is target.fa.
### 1. Genome alignment preprocessing
```bash
./01_pre.sh -s species.txt -r genome1 -d /Yuli-Lab-wangyajun/DistantSV-Cross/example

*Parameter description for 01_pre.sh*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference
   -d working_directory  : Working directory path
```
Note: This step will automatically create four directories in the specified working path "/Yuli-Lab-wangyajun/DistantSV-Cross/example": ref_chr, which stores the reference genome split by sequence count to enhance alignment efficiency; species_2bit, which generates 2bit files for each species genome for subsequent multiple sequence alignment; species_nibs, which produces nib files for each species genome; and species_sizes, which generates sizes files for each species genome. Specifically, it is recommended to use a reference genome assembled to the chromosome level, although a scaffold-level genome is also acceptable; however, a higher number of sequences in the reference genome will increase the task count, though this does not affect the analysis. Additionally, if there are many unassembled short sequence files, you can manually delete unnecessary ones in the ref_chr directory, retaining only the sequence files required for alignment.

### 2. Pairwise genome alignment
```bash
./02_align.sh -s species.txt -r genome1 -d /Yuli-Lab-wangyajun/DistantSV-Cross/example -t Template1_lastz_script.sh

*Parameter description for 02_align.sh*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference
   -d working_directory  : Working directory path
   -t template.sh        : LASTZ template script
```
```bash
*Template1_lastz_script.sh*
#!/bin/bash
#SBATCH -J JOB_NAME_PLACEHOLDER
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o /work/home/acm2g39zj0/err_out/std.out.%j
#SBATCH -e /work/home/acm2g39zj0/err_out/std.err.%j

source ~/miniconda3/bin/activate
conda activate env_test
cd WORK_DIR_PLACEHOLDER/axt/CHROMOSOME_PLACEHOLDER_axt
lastz WORK_DIR_PLACEHOLDER/ref_chr/CHROMOSOME_PLACEHOLDER GENOME_NAME_PLACEHOLDER --notransition --ydrop=15000 --gappedthresh=3000 --hspthresh=4500 --format=axt --out=GENOME_NAME_OUTPUT.axt
```
Note: This step involves batch-generating alignment scripts. To run the commands mentioned above, you need to modify the "#SBATCH" parameters at the beginning of the lastz_template.sh template file provided by us, based on the job scheduling system of your cluster. The lastz alignment parameters, such as "--notransition --ydrop=15000 --gappedthresh=3000 --hspthresh=4500 --format=axt," can be adjusted as needed. However, this set of parameters has been evaluated by us as the most suitable for aligning distantly related genomes. Therefore, we recommend using the default parameters in the template file. No other modifications are required elsewhere in the file.This step will automatically generate two directories in the specified working path "/Yuli-Lab-wangyajun/DistantSV-Cross/example":
shell_sh: Contains batch-generated alignment scripts for sequentially aligning each other genome to each chromosome of the reference genome; axt: Stores the pairwise alignment results for each chromosome of the reference genome against the other genomes. To facilitate efficient batch submission, we have generated a file named batch_lastz_submit.sh, which contains all submission scripts. You can simply run “./batch_lastz_submit.sh” to submit all tasks at once. However, we recommend submitting tasks in batches from the shell_sh directory based on the actual resource availability of your cluster.

### 3. Generating Progressive Multi-Sequence Alignment Files
```bash
./MSA.sh -s species.txt -r genome1 -d /Yuli-Lab-wangyajun/DistantSV-Cross/example -t Template2_multiz_script.sh -tree "(((genome1 genome2) target) genome3)"

*Parameter description for 03_MSA.sh*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference
   -d working_directory  : Working directory path
   -t template.sh        : MULTIZ template script
   -tree "species_tree"  : species tree
```
```bash
*Template2_multiz_script.sh*
#!/bin/bash 
#SBATCH -J JOB_NAME_PLACEHOLDER
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o /work/home/acm2g39zj0/err_out/multiz_CHROMOSOME_PLACEHOLDER.out
#SBATCH -e /work/home/acm2g39zj0/err_out/multiz_CHROMOSOME_PLACEHOLDER.err

source ~/miniconda3/bin/activate
conda activate env_test
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
```
Note: This step involves batch-generating multiple sequence alignment scripts. You need to modify the "#SBATCH" parameters at the beginning of the Template2_multiz_script.sh template file provided by us, based on the job scheduling system of your cluster. No other modifications are required elsewhere in the file.This step will batch-generate scripts for constructing multiple sequence alignments for each chromosome of the reference genome in the shell_sh directory. Similarly, users need to submit tasks in batches from the shell_sh directory based on the actual resource availability of their cluster.Before submitting the tasks, the pairwise alignments between species in step2 must be completed.After the multiple sequence alignment tasks are completed in batches, a directory for multiple sequence alignment results will be automatically created in the working directory "/Yuli-Lab-wangyajun/DistantSV-Cross/example". The multiz_out directory will contain multiple sequence alignment files with the .maf suffix. This step completes the multiple sequence alignment of each sequence file of the reference genome with the genomes of other species. Subsequent analysis will process the MAF-format multiple sequence alignment files and ultimately detect structural variations in the target species.


### 4-5-6. Alignment block filtering and realignment
```bash
./SV.sh -s species.txt -r genome1 -d /Yuli-Lab-wangyajun/DistantSV-Cross/example -t Template3_SV_script.sh -target target -non_target "genome1,genome2,genome3"

*Parameter description for 04_SV.sh*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference species
   -d working_directory  : Working directory path
   -t template.sh        : SV template script
   -target target_species: Target species
   -non_target "list"    : Non-target species list,‘,’delimited
```
```bash
*Template3_SV_script.sh*
#!/bin/bash
#SBATCH -J JOB_NAME_PLACEHOLDER
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH -o /work/home/acm2g39zj0/err_out/std.out.%j
#SBATCH -e /work/home/acm2g39zj0/err_out/std.err.%j

source ~/miniconda3/bin/activate
conda activate env_test
cd WORK_DIR_PLACEHOLDER/CHROMOSOME_PLACEHOLDER_multiz_out/multiz_out
python ../../Module1_MAF_filter.py --input CHROMOSOME_PLACEHOLDER.maf --output CHROMOSOME_PLACEHOLDER_filter.maf --min_sequences 4 --ref_species REF_GENOME_PLACEHOLDER --min_ref_length 500 --max_gap_ratio 0.1 --gap_region_length 100
python ../../Module2_MAF_realigner.py -i CHROMOSOME_PLACEHOLDER_filter.maf -o CHROMOSOME_PLACEHOLDER_filter_realign.maf
python ../../Module3_call_SV.py --maf CHROMOSOME_PLACEHOLDER_filter_realign.maf --target TARGET_SPECIES_PLACEHOLDER --reference REF_GENOME_PLACEHOLDER --non_target "NON_TARGET_PLACEHOLDER" --output_prefix CHROMOSOME_PLACEHOLDER --max_target_gap 10 --min_ref_gap 50 --min_target_gap 50 --max_ref_gap 10 --max_gap 100 --inv_more_block_max 1000 --min_intra_size 50 --max_gap_ratio 0.3
python ../../Module4_sequence_consistency_calculation.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_insertions.bed -o CHROMOSOME_PLACEHOLDER_ins.bed -r REF_GENOME_PLACEHOLDER
python ../../Module4_sequence_consistency_calculation.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_deletions.bed -o CHROMOSOME_PLACEHOLDER_del.bed -r REF_GENOME_PLACEHOLDER
python ../../Module4_sequence_consistency_calculation.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_inversions.bed -o CHROMOSOME_PLACEHOLDER_inv.bed -r REF_GENOME_PLACEHOLDER
python ../../Module4_sequence_consistency_calculation.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_interchrom_translocations.bed -o CHROMOSOME_PLACEHOLDER_interTRA.bed -r REF_GENOME_PLACEHOLDER
python ../../Module4_sequence_consistency_calculation.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_intrachrom_translocations.bed -o CHROMOSOME_PLACEHOLDER_intraTRA.bed -r REF_GENOME_PLACEHOLDER
```
Note: This step integrates three steps—"step4: Alignment block filtering and realignment," "step5: Calling structural variation," and "step6: Calculation of the Mean Sequence Consistency"—into one, batch-generating SV detection scripts. Meanwhile, you need to modify the Template3_SV_script.sh template file provided by us, specifically adjusting the "#SBATCH" parameters at the beginning based on the job scheduling system of your cluster. Additionally, in the provided template file Template3_SV_script.sh, you must modify some parameters in "step4: Calling structural variation," such as "--min_sequences 4" (the example includes 4 species; you need to adjust this based on the number of species in your dataset). For other steps, it is recommended to keep the default parameters. Finally, we will also provide detailed parameter usage for [Module1_MAF_filter.py](#Module1_MAF_filter.py), [Module2_MAF_realigner.py](#Module2_MAF_realigner.py), [Module3_call_SV.py](#Module3_call_SV.py), and [Module4_sequence_consistency_calculation.py](#Module4_sequence_consistency_calculation.py) below, allowing you to make flexible adjustments. After running the command, scripts will be batch-generated in the shell_sh directory for each multiple sequence alignment MAF file. You can submit jobs as described earlier. The final results will be generated in the multiz_out directory where each MAF file is located, producing SV files in [BED format](#BED). The final result files are BED files ending with "_ins.bed", "_del.bed", "_inv.bed", "_intraTRA.bed", and "_interTRA.bed". Users only need to focus on these five result files.

## Detailed parameters

Parameter description for <a id="Module1_MAF_filter.py"></a>Module1_MAF_filter.py
```
Function: Filters multiple sequence alignment blocks constructed by MULTIZ to exclude inappropriate alignment blocks

Usage: python Module1_MAF_filter.py [-h] --input INPUT --output OUTPUT [--min_sequences MIN_SEQUENCES] [--required_species REQUIRED_SPECIES] [--ref_species REF_SPECIES] [--min_ref_length MIN_REF_LENGTH] [--max_gap_ratio MAX_GAP_RATIO] [--gap_region_length GAP_REGION_LENGTH] [--species_gap_ratio SPECIES_GAP_RATIO] [--species_gap_region_length SPECIES_GAP_REGION_LENGTH]

Options:
-h, --help show this help message and exit
--input INPUT (required) Path to the input MAF file.
--output OUTPUT (required) Path to the output MAF file.
--min_sequences MIN_SEQUENCES
Minimum number of sequences (species) required in a block. (default: 2)
--required_species REQUIRED_SPECIES
Comma-separated list of species that must appear in every block. (optional)
--ref_species REF_SPECIES
Reference species name used for length checks. (default: ho_sa)
--min_ref_length MIN_REF_LENGTH
Minimum alignment length of the reference sequence. (default: 100)
--max_gap_ratio MAX_GAP_RATIO
Maximum allowed gap ratio at the left and right ends of the alignment (global check). (default: 0.05)
--gap_region_length GAP_REGION_LENGTH
Length of the end regions (in bases) considered for global gap check. (default: 100)
--species_gap_ratio SPECIES_GAP_RATIO
Species-specific gap thresholds, e.g., "species1:0.1,species2:0.2". (optional)
--species_gap_region_length SPECIES_GAP_REGION_LENGTH
Length of end regions for species-specific gap checks. If not set, uses --gap_region_length. (optional)
```

Parameter description for <a id=Module2_MAF_realigner.py></a>Module2_MAF_realigner.py
```
Function: Performs multiple sequence realignment on the filtered multiple sequence alignment blocks using MUSCLE and MAFFT.

Usage: python Module2_MAF_realigner.py [-h] -i INPUT -o OUTPUT [--muscle MUSCLE] [--mafft MAFFT] [-t THRESHOLD] [-v]

Options:
-h, --help show this help message and exit
-i INPUT, --input INPUT
(required) Input MAF file.
-o OUTPUT, --output OUTPUT
(required) Output MAF file.
--muscle MUSCLE Path to MUSCLE5 executable. (default: muscle)
--mafft MAFFT Path to MAFFT executable. (default: mafft)
-t THRESHOLD, --threshold THRESHOLD
Sequence length threshold: use MAFFT if above, otherwise MUSCLE5. (default: 1000)
-v, --verbose Verbose output (set logging to DEBUG).
```

Parameter description for <a id=Module3_call_SV.py></a>Module3_call_SV.py
```
Function: Perform SV detection on the realigned MAF files.

Usage: python Module3_SV_detection.py [-h] --maf MAF --target TARGET --reference REFERENCE --output_prefix OUTPUT_PREFIX [--non_target NON_TARGET] [--max_target_gap MAX_TARGET_GAP] [--min_ref_gap MIN_REF_GAP] [--min_target_gap MIN_TARGET_GAP] [--max_ref_gap MAX_REF_GAP] [--max_gap MAX_GAP] [--inv_more_block_max INV_MORE_BLOCK_MAX] [--min_intra_size MIN_INTRA_SIZE] [--max_gap_ratio MAX_GAP_RATIO]

Options:
-h, --help show this help message and exit
--maf MAF (required) Input MAF file.
--target TARGET (required) Target species name.
--reference REFERENCE (required) Reference species name.
--output_prefix OUTPUT_PREFIX
(required) Output file prefix.
--non_target NON_TARGET
Other non-target species (comma-separated). (default: "")
--max_target_gap MAX_TARGET_GAP
Maximum continuous gap for target species, used in deletion detection. (default: 10)
--min_ref_gap MIN_REF_GAP
Minimum SV gap for reference species. (default: 50)
--min_target_gap MIN_TARGET_GAP
Minimum SV gap for target species. (default: 50)
--max_ref_gap MAX_REF_GAP
Maximum continuous gap for reference species, used in insertion detection. (default: 10)
--max_gap MAX_GAP Maximum gap for continuous regions, used in translocation and inversion detection. (default: 100)
--inv_more_block_max INV_MORE_BLOCK_MAX
Maximum allowed gap threshold for all species in multi‑block inversion detection. (default: 1000)
--min_intra_size MIN_INTRA_SIZE
Minimum size for intra‑block SVs. (default: 10)
--max_gap_ratio MAX_GAP_RATIO
Maximum allowed gap ratio for non‑target species in deletion detection. (default: 0.2)
```

Parameter description for <a id=Module4_sequence_consistency_calculation.py></a>Module4_sequence_consistency_calculation.py
```
Function: Calculate the mean sequence consistency of the flanking sequences on both sides of the SV.

Usage: Module4_conservation_calc.py [-h] -m MAF -b BED -o OUTPUT [-r REF_SPECIES]

Options:
-h, --help show this help message and exit
-m MAF, --maf MAF (required) Input MAF file.
-b BED, --bed BED (required) Input BED file.
-o OUTPUT, --output OUTPUT
(required) Output BED file.
-r REF_SPECIES, --ref-species REF_SPECIES
Reference species name (default: mu_mu)
```

<a id=BED></a>BED format
```
We have generated a total of 22 columns of variation information, providing detailed records of the detected target species-specific structural variations in both the target genome and the corresponding reference genome coordinates, offering users diverse criteria for screening candidate SVs.

Ref_Chrom: Reference genome chromosome ID
Ref_Start: Reference genome start position
Ref_End: Reference genome end position
Ref_Strand: Reference genome strand orientation information
Target_Chrom: Target genome chromosome ID
Target_Start: Target genome start position
Target_End: Target genome end position
Target_Strand: Target genome strand orientation information
Original_tra_Chrom: Chromosome ID of the translocated fragment, only present for translocation variations
Original_tra_Start: Genome start position of the translocated fragment, only present for translocation variations
Original_tra_End: Genome end position of the translocated fragment, only present for translocation variations
Original_tra_Strand: Strand orientation information of the translocated fragment, only present for translocation variations
SV_Type: Type of genomic structural variation
SV_Size: Size of the genomic structural variation (bp)
SV_Source: Context of the variation within the alignment block, such as whether it is inside or between alignment blocks
Block_Indexes: Indexes of alignment blocks involved in the variation
Additional_Info: Additional information, same as SV_Source
SV_Upstream_Sequence_Consistency: Sequence consistency upstream of the variation position
SV_Downstream_Sequence_Consistency: Sequence consistency downstream of the variation position
SV_Total_Sequence_Consistency: Total sequence consistency upstream and downstream of the variation position
SV_Upstream_Length: Calculation range for sequence consistency upstream of the variation position
SV_Downstream_Length: Calculation range for sequence consistency downstream of the variation position
```

## License
This project is protected by the GNU General Public License v3.0 license.
