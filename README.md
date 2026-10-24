## Overview
<div align="center"><img width="85%" src="https://github.com/Yuli-Lab-wangyajun/DistantSV-Cross/blob/master/picture/DistantSV-Cross-logo.png"/></div>

DistantSV-Cross is a computational pipeline specifically designed to detect structural variations (SVs) among distantly related multiple species. It consists of six computational steps to identify species-specific genomic structural variations in a target species relative to others, ultimately generating information on insertions, deletions, inversions, intra-chromosomal translocations, and inter-chromosomal translocations in BED format. Compared to other variant detection tools based on pairwise species comparisons via genome assembly, as well as existing multi-species structural variation detection software, DistantSV-Cross demonstrates superior performance, achieving higher precision, recall, and F1 scores. <br/>

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
  - [3. Progressive multi-sequence alignment](#3-progressive-multi-sequence-alignment)
  - [4-5-6. Alignment block filtering and realignment; Calling SVs; Mean sequence consistency of SV breakpoints](#4-5-6-alignment-block-filtering-and-realignment-calling-svs-mean-sequence-consistency-of-sv-breakpoints)
  - [Serial and parallel computing strategies](#serial-and-parallel-computing-strategies)
- [Detailed parameters](#detailed-parameters)
- [License](#license)
- [Citing DistantSV-Cross](#citing-distantsv-cross)
- [Gitting help](#gitting-help)

## Installation:
DistantSV-Cross offers a comprehensive pipeline for running the code, requiring only the installation of the necessary software, Python packages, and downloading the source code provided by us.

### 1. Download the source code
```
git clone https://github.com/Yuli-Lab-wangyajun/DistantSV-Cross
or
If you encounter "Connection timed out", you can also manually download the source code from the "Code" directory.
```

### 2. Install the required dependencies for running the software
```
## We recommend creating a new environment (e.g., test) with Python 3.8 or higher to install the required software or packages.
conda create -n test python=3.8 pandas numpy -y

## Install the following software manually and make sure they are accessible within the "test" environment.
lastz, multiz, muscle, mafft

## Download the following packages from https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ and add them to the env_test environment.
axtChain, chainPreNet, chainNet, netToAxt, axtToMaf, faSize, faToTwoBit, faSplit
```

## Data preparation
To demonstrate the usage of DistantSV-Cross, we will use the sample files provided in the demo folder. All steps described below should be executed in the working directory “/Yuli-Lab-wangyajun/DistantSV-Cross/example_supplementary_data/example”, which is contained in the example_supplementary_data.zip file in the released “V1.0”.

Before running, the following information needs to be prepared:
### 1. Prepare the species.txt file
```
genome1	/Yuli-Lab-wangyajun/DistantSV-Cross/example_supplementary_data/example/genome/genome1.fa
genome2	/Yuli-Lab-wangyajun/DistantSV-Cross/example_supplementary_data/example/genome/genome2.fa
genome3	/Yuli-Lab-wangyajun/DistantSV-Cross/example_supplementary_data/example/genome/genome3.fa
target	/Yuli-Lab-wangyajun/DistantSV-Cross/example_supplementary_data/example/genome/target.fa
```
Note: The first column specifies an abbreviated name for each species' genome, which you can define according to your needs—keeping it concise is recommended. The second column must contain the absolute path to the corresponding genome file. Please ensure that these two columns are separated by a tab. In this example, target.fa is designated as the target species to be analyzed, while genome1.fa serves as the reference genome.

### 2. Determine the phylogenetic relationships among the species
```
(((genome1 target) genome2) genome3);
```
Note: The phylogenetic relationships in Newick format provide the input for constructing multiple sequence alignments with MULTIZ in subsequent steps.

## Usage 
With the above information prepared, we can now initiate the SV detection process. A fully functional working example is provided in the "example" directory, and the step-by-step operational procedure is detailed below:
<div align="center"><img width="95%" src="https://github.com/Yuli-Lab-wangyajun/DistantSV-Cross/blob/master/picture/graphical_abstract.png"/></div>
### 1. Genome alignment preprocessing
```bash
./01_pre.sh -s species.txt -r genome1 -d /work/home/acm2g39zj0/wangyj/SV/new_SV_method/example

*Parameter description for 01_pre.sh*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference
   -d working_directory  : Working directory path
```
Note: This step automatically creates four directories under the specified working path: (1) ref_chr: stores the reference genome split by sequence count to improve alignment efficiency; (2) species_2bit: contains 2bit files for each species genome, used for subsequent multiple sequence alignment; (3) species_nibs: holds nib files for each species genome; (4) species_sizes: stores size files for each species genome.<br/>
It is recommended to use a chromosome-level reference genome for optimal performance, although a scaffold-level genome is acceptable. Note that a reference genome with a larger number of sequences will increase the number of tasks, but this does not affect the analysis results. If the ref_chr directory contains many unassembled short sequence files, you may manually delete unnecessary ones and retain only the sequences required for alignment.

### 2. Pairwise genome alignment
```bash
./02_align.sh -s species.txt -r genome1 -d /Yuli-Lab-wangyajun/DistantSV-Cross/example_supplementary_data/example -t Template1_lastz_script.sh

*Parameter description for 02_align.sh*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference
   -d working_directory  : Working directory path
   -t template.sh        : LASTZ template script
```
Note: This step involves batch-generating LASTZ alignment scripts. To run the commands described above, you need to modify the #SBATCH parameters at the beginning of the Template1_lastz_script.sh template file (provided by us) according to the job scheduling system of your cluster.This step will automatically create two directories under the specified working path: (1) shell_sh: Contains the batch-generated alignment scripts for sequentially aligning each other genome to each chromosome of the reference genome; (2) axt: Stores the pairwise alignment results for each chromosome of the reference genome against the other genomes.<br/>
To facilitate efficient batch submission, we have also generated a file named batch_lastz_submit.sh, which contains all submission commands.You can run "./batch_lastz_submit.sh" to submit all tasks at once, or choose to submit the parallel task "combined_XXX_job.sh". However, we recommend submitting tasks in batches from the shell_sh directory, depending on the actual resource availability of your cluster.

### 3. Progressive multi-sequence alignment
```bash
./03_MSA.sh -s species.txt -r genome1 -d /Yuli-Lab-wangyajun/DistantSV-Cross/example_supplementary_data/example -t Template2_multiz_script.sh -tree "(((genome1 genome2) target) genome3)"

*Parameter description for 03_MSA.sh*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference
   -d working_directory  : Working directory path
   -t template.sh        : MULTIZ template script
   -tree "species_tree"  : species tree
```
Note: This step involves batch-generating multiple sequence alignment scripts. You need to modify the "#SBATCH" parameters at the beginning of the Template2_multiz_script.sh template file provided by us, based on the job scheduling system of your cluster. No other modifications are required elsewhere in the file.This step will batch-generate scripts for constructing multiple sequence alignments for each chromosome of the reference genome in the shell_sh directory. Similarly, users need to submit tasks in batches from the shell_sh directory based on the actual resource availability of their cluster.Before submitting the tasks, the pairwise alignments between species in step2 must be completed.After the multiple sequence alignment tasks are completed in batches, a directory for multiple sequence alignment results will be automatically created in the working directory. The multiz_out directory will contain multiple sequence alignment files with the .maf suffix. This step completes the multiple sequence alignment of each sequence file of the reference genome with the genomes of other species. Subsequent analysis will process the MAF-format multiple sequence alignment files and ultimately detect structural variations in the target species.<br/>
This step involves batch-generating multiple sequence alignment scripts. You need to modify the #SBATCH parameters at the beginning of the Template2_multiz_script.sh template file according to the job scheduling system of your cluster. No other modifications to the file are required.<br/>
Before submitting the alignment tasks, ensure that the pairwise alignments between species from "2. Pairwise genome alignment" have been successfully completed.<br/>
This step will automatically generate scripts for constructing multiple sequence alignments for each chromosome of the reference genome, which will be saved in the shell_sh directory. You can run "./batch_multiz_submit.sh" to submit all tasks at once, or choose to submit the parallel task "combined_multiz_job.sh". However, we recommend submitting tasks in batches from the shell_sh directory, depending on the actual resource availability of your cluster.<br/>
Once the multiple sequence alignment tasks are completed, multiple directories ending with "_multiz_out" will be automatically generated under the working path, each of which contains multiple sequence alignment result files (with the .maf suffix).

### 4-5-6. Alignment block filtering and realignment; Calling SVs; Mean sequence consistency of SV breakpoints
```bash
./04_SV.sh -s species.txt -r genome1 -d /Yuli-Lab-wangyajun/DistantSV-Cross/example_supplementary_data/example -t Template3_SV_script.sh -target target -non_target "genome1,genome2,genome3"

*Parameter description for 04_SV.sh*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference species
   -d working_directory  : Working directory path
   -t template.sh        : SV template script
   -target target_species: Target species
   -non_target "list"    : Non-target species list,‘,’delimited
```
Note: This step integrates three steps "4: Alignment block filtering and realignment," "5: Calling SVs," and "6: Mean sequence consistency of SV breakpoints" into one, batch-generating SV detection scripts. Meanwhile, you need to modify the Template3_SV_script.sh template file provided by us, specifically adjusting the "#SBATCH" parameters at the beginning based on the job scheduling system of your cluster. Additionally, in the provided template file Template3_SV_script.sh, you must modify some parameters in "step4: Calling structural variation," such as "--min_sequences 4" (the example includes 4 species; you need to adjust this based on the number of species in your dataset). For other steps, it is recommended to keep the default parameters. <br/>
Finally, we will also provide detailed parameter usage for [Module1_MAF_filter.py](#Module1_MAF_filter.py), [Module2_MAF_realigner.py](#Module2_MAF_realigner.py), [Module3_call_SV.py](#Module3_call_SV.py), and [Module4_sequence_consistency_calculation.py](#Module4_sequence_consistency_calculation.py) below, allowing you to make flexible adjustments. After running the command, scripts will be batch-generated in the shell_sh directory for each multiple sequence alignment MAF file. You can run "./batch_SV_submit.sh" to submit all tasks at once, or choose to submit the parallel task "combined_SV_job.sh". However, we recommend submitting tasks in batches from the shell_sh directory, depending on the actual resource availability of your cluster. The final results will be generated in the multiz_out directory where each MAF file is located, producing SV files in [BED format](#BED). The final result files are BED files ending with "_ins.bed", "_del.bed", "_inv.bed", "_intraTRA.bed", and "_interTRA.bed". Users only need to focus on these five result files.<br/>
The example folder contains a simulated file named Simulated_SV.bed, which includes 10 SVs (2 INS, 2 DEL, 2 INV, 2 Intra-TRA, and 2 Inter-TRA). The actual detection results are also provided in the example folder for review.

### Serial and parallel computing strategies
```bash
python GNUparallel.py /Yuli-Lab-wangyajun/DistantSV-Cross/example_supplementary_data/example/shell_sh -s genome1,genome2,genome3,target -r genome1 -g 60
```
Note: Based on the GNU Parallel tool, we additionally provide parallel task execution. Users can directly submit combined_*_job.sh to achieve parallel execution of a single task, or choose to run multiple tasks simultaneously, and select the appropriate running mode based on available computing resources.<br/>

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

## Citing DistantSV-Cross
If you use DistantSV-Cross in your research, please cite our repository (and forthcoming paper):<br/>

(https://github.com/Yuli-Lab-wangyajun/DistantSV-Cross)
<br/>

## Gitting help
If you encounter any issues or have questions about specific parameters, please Open an Issue or contact with Yajun Wang (yajunwang2024@foxmail.com).
<br/>
