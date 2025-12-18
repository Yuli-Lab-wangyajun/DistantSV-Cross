DistantSV-Cross is a computational pipeline specifically designed for detecting structural variations among distantly related multiple species. It primarily comprises six computational steps to identify species-specific genomic structural variations in a target species relative to other species, ultimately generating insertion, deletion, inversion, intra-chromosomal translocation, and inter-chromosomal translocation information in BED format. Compared to other variant detection tools based on pairwise species comparison via genome assembly, as well as existing multi-species structural variation detection software, DistantSV-Cross demonstrates superior performance, achieving higher precision, recall, and F1-score. <br/>

## Installation:
```
## 建议新建一个环境安装所依赖的软件或包，python环境在3.8+，这里以env_test为例
conda create -n env_test python=3.8 pandas numpy -y

## 手动安装下列软件，并添加到env_test环境中
这四个可单独安装并手动添加到env_test环境：
lastz、multiz、muscle、mafft
来自（https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/）的安装包，寻找下列安装包并添加到env_test环境：
axtChain、chainPreNet、chainNet、netToAxt、axtToMaf、faSize、faToTwoBit、faSplit
```

## 需要提前准备的信息
我这里以example文件夹中的示例文件来详细说明DistantSV-Cross的使用步骤，下述所有的操作步骤均应该在“/Yuli-Lab-wangyajun/DistantSV-Cross/example”工作路径下进行
### 准备1: species.txt文件的准备
```
*species.txt*
genome1	/Yuli-Lab-wangyajun/DistantSV-Cross/example/genomes_path/genome1.fa
genome2	/Yuli-Lab-wangyajun/DistantSV-Cross/example/genomes_path/genome2.fa
genome3	/Yuli-Lab-wangyajun/DistantSV-Cross/example/genomes_path/genome3.fa
target	/Yuli-Lab-wangyajun/DistantSV-Cross/example/genomes_path/target.fa
```
第一列是每个物种基因组的命名简称，方便你自己辨识即可，建议简短一些。第二列是基因组的绝对路径。注意，这两列之间要用tab键分割。这个示例中是确定target.fa为我们待检测的目标物种，参考物种为genome1.fa。

### 准备2: 确定好物种之间的系统发育关系
```
*nwk格式的系统发育关系*
(((genome1 target) genome2) genome3);
```
Newick格式的系统发育关系，主要用于后续multiz构建多序列比对文件时使用

### 准备3: 拷贝下列工作脚本到工作目录
```
01_align.sh
maf_filter.py
02_MSA.sh
new_bao.py
species.txt
sv_template.sh
all.sh
conservation_calculator.py
lastz_template.sh
MAF_realigner.py
multiz_template.sh
00_pre.sh
SV.sh
```

## 具体运行步骤:
### step1: Genome alignment preprocessing:
```
*命令行*
./pre.sh -s species.txt -r genome1 -d /Yuli-Lab-wangyajun/DistantSV-Cross/example

*pre.sh的参数描述*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference
   -d working_directory  : Working directory path
```

此步骤将在指定的工作路径“/Yuli-Lab-wangyajun/DistantSV-Cross/example”下自动生成4个目录：ref_chr（将参考基因组按照序列条数进行拆分，目的是为提高后续比对效率）、species_2bit（为每个物种基因组生成2bit文件，以便后续构建多序列比对文件）、species_nibs（为每个物种基因组生成nib文件，以便后续构建多序列比对文件）、species_sizes（为每个物种基因组生成sizes文件，以便后续构建多序列比对文件）。需要注意的是，建议选择一个挂载到染色体水平的参考基因组，也可以是scaffold水平，只是当参考基因组中序列条数变多之后会使得运行任务数量的增多，不过这并不影响分析。如果基因组未挂载的短序列文件较多，你可手动删掉这些序列文件，保留较长的序列文件。

### spet2: Pairwise genome alignment:
```
*命令行*
./align.sh -s species.txt -r genome1 -d /Yuli-Lab-wangyajun/DistantSV-Cross/example -t lastz_template.sh

*align.sh的参数描述*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference
   -d working_directory  : Working directory path
   -t template.sh        : LASTZ template script
```
```
*lastz_template.sh*
#!/bin/bash
#SBATCH -J JOB_NAME_PLACEHOLDER
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o /work/home/acm2g39zj0/err_out/std.out.%j
#SBATCH -e /work/home/acm2g39zj0/err_out/std.err.%j

conda activate env_test
cd WORK_DIR_PLACEHOLDER/axt/CHROMOSOME_PLACEHOLDER_axt
lastz WORK_DIR_PLACEHOLDER/ref_chr/CHROMOSOME_PLACEHOLDER GENOME_NAME_PLACEHOLDER --notransition --ydrop=15000 --gappedthresh=3000 --hspthresh=4500 --format=axt --out=GENOME_NAME_OUTPUT.axt
```
想要运行上述命令，你还需要在我们提供的lastz_template.sh模板文件，根据你所在集群作业调度系统，修改开头的“#SBATCH”参数部分即可。lastz比对参数“--notransition --ydrop=15000 --gappedthresh=3000 --hspthresh=4500 --format=axt”可自行修改，我们建议用模板文件中的默认参数即可，其余地方均不用修改。此步骤将在指定的工作路径“/Yuli-Lab-wangyajun/DistantSV-Cross/example”下自动生成2个目录：shell_sh（批量生成其他基因组依次比对到参考基因组的每一条染色体上的比对脚本）、axt（用于存放参考基因组的每一条染色体的两两物种比对结果）。step2完成后，，为方便快速批量提交，我们生成了提交所有脚本的文件batch_lastz_submit.sh，你只需运行“./batch_lastz_submit.sh”即可批量全部提交。需要注意的是，我们建议根据你自己实际的集群资源情况在shell_sh目录下分批次提交任务。

### step3: Generating Progressive Multi-Sequence Alignment Files:
```
*命令行*
./MSA.sh -s species.txt -r genome1 -d /Yuli-Lab-wangyajun/DistantSV-Cross/example -t multiz_template.sh -tree "(((genome1 genome2) target) genome3)"

*MSA.sh的参数描述*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference
   -d working_directory  : Working directory path
   -t template.sh        : MULTIZ template script
   -tree "species_tree"  : species tree (required)
```
```
*multiz_template.sh*
#!/bin/bash 
#SBATCH -J JOB_NAME_PLACEHOLDER
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o /work/home/acm2g39zj0/err_out/multiz_CHROMOSOME_PLACEHOLDER.out
#SBATCH -e /work/home/acm2g39zj0/err_out/multiz_CHROMOSOME_PLACEHOLDER.err

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
与step2一致，此步骤也需要在我们提供的multiz_template.sh模板文件，根据你所在集群作业调度系统，修改开头的“#SBATCH”参数部分即可，其他地方均不用修改。此步骤将在shell_sh目录下为参考基因组的每一条染色体批量生成多序列比对构建的脚本。同理，用户也需要根据自己实际的集群资源情况在shell_sh目录下分批次提交任务。在提交任务之前，必须先完成step2中的两两物种间的比对。多序列比对任务分批次运行结束后，将在工作目录“/Yuli-Lab-wangyajun/DistantSV-Cross/example”下自动建立并生成多序列比对结果目录，在multiz_out目录中存放着以maf为后缀的多序列比对文件，这一步便完成了参考基因组的每一条序列文件与其他物种基因组的多序列比对结果，后续分析将对MAF格式的多序列比对文件进行处理以及最终目标物种结构变异的检测。

### step4-5-6: Alignment block filtering and realignment
```
*命令行*
./SV.sh -s species.txt -r genome1 -d /Yuli-Lab-wangyajun/DistantSV-Cross/example -t sv_template.sh -target target -non_target "genome1,genome2,genome3"

*SV.sh的参数描述*
   -s species.txt        : Files of genome abbreviations and their paths
   -r reference_genome   : Reference species
   -d working_directory  : Working directory path
   -t template.sh        : SV template script
   -target target_species: Target species
   -non_target "list"    : Non-target species list,‘,’delimited
```
```
*sv_template.sh*
#!/bin/bash
#SBATCH -J JOB_NAME_PLACEHOLDER
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o /work/home/acm2g39zj0/err_out/std.out.%j
#SBATCH -e /work/home/acm2g39zj0/err_out/std.err.%j

conda activate env_test
cd WORK_DIR_PLACEHOLDER/CHROMOSOME_PLACEHOLDER_multiz_out/multiz_out
##step4: Alignment block filtering and realignment
python ../../maf_filter.py --input CHROMOSOME_PLACEHOLDER.maf --output CHROMOSOME_PLACEHOLDER_filter.maf --min_sequences 4 --ref_species REF_GENOME_PLACEHOLDER --min_ref_length 100 --max_gap_ratio 0.3 --gap_region_length 50
python ../../MAF_realigner.py -i CHROMOSOME_PLACEHOLDER_filter.maf -o CHROMOSOME_PLACEHOLDER_filter_realign.maf
##step5: Calling structural variation
python ../../sv_Calling.py --maf CHROMOSOME_PLACEHOLDER_filter_realign.maf --target TARGET_SPECIES_PLACEHOLDER --reference REF_GENOME_PLACEHOLDER --non_target "NON_TARGET_PLACEHOLDER" --output_prefix CHROMOSOME_PLACEHOLDER --max_target_gap 10 --min_ref_gap 50 --min_target_gap 50 --max_ref_gap 10 --max_gap 100 --inv_more_block_max 1000 --min_intra_size 50 --max_gap_ratio 0.3
##step6: Calculation of the Mean Sequence Consistency
python ../../conservation_calculator.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_insertions.bed -o CHROMOSOME_PLACEHOLDER_ins.bed -r REF_GENOME_PLACEHOLDER
python ../../conservation_calculator.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_deletions.bed -o CHROMOSOME_PLACEHOLDER_del.bed -r REF_GENOME_PLACEHOLDER
python ../../conservation_calculator.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_inversions.bed -o CHROMOSOME_PLACEHOLDER_inv.bed -r REF_GENOME_PLACEHOLDER
python ../../conservation_calculator.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_interchrom_translocations.bed -o CHROMOSOME_PLACEHOLDER_interTRA.bed -r REF_GENOME_PLACEHOLDER
python ../../conservation_calculator.py -m CHROMOSOME_PLACEHOLDER_filter_realign.maf -b CHROMOSOME_PLACEHOLDER_intrachrom_translocations.bed -o CHROMOSOME_PLACEHOLDER_intraTRA.bed -r REF_GENOME_PLACEHOLDER
```
与step2一致，此步骤也需要在我们提供的sv_template.sh模板文件，根据你所在集群作业调度系统，修改开头的“#SBATCH”参数部分即可。命令运行结束后将在shell_sh目录下为每一个多序列比对MAF文件批量生成为生成脚本。脚本将3个步骤“step4: Alignment block filtering and realignment; step5: Calling structural variation; step6: Calculation of the Mean Sequence Consistency”整合在一起了。需要注意的是，你需要根据我们提供的模板文件sv_template.sh中，切记要修改“step4: Calling structural variation”中的部分参数，如“--min_sequences 4”，示例中有4个物种，你需要根据自己数据集的物种数目进行修改，其他步骤建议默认参数即可，我们也将在下面将分别补充[maf_filter.py](#maf_filter.py)、[MAF_realigner.py](#MAF_realigner.py)、[sv_Calling.py](#sv_Calling.py)、[conservation_calculator.py](#conservation_calculator.py)参数的详细情况。结果会在每一个多序列MAF文件所在的multiz_out目录下生成[bed格式的SV文件](#bed格式的SV文件)，最终结果文件是以“_ins.bed”、“_del.bed”、“_inv.bed”、“_intraTRA.bed”、“_interTRA.bed”结尾的bed文件，用户只需关注这五个结果文件即可。

## 补充说明

<a name=maf_filter.py></a>maf_filter.py的参数说明
```
usage: maf_filter.py [-h] --input INPUT --output OUTPUT [--min_sequences MIN_SEQUENCES] [--required_species REQUIRED_SPECIES] [--ref_species REF_SPECIES] [--min_ref_length MIN_REF_LENGTH]
                     [--max_gap_ratio MAX_GAP_RATIO] [--gap_region_length GAP_REGION_LENGTH] [--species_gap_ratio SPECIES_GAP_RATIO] [--species_gap_region_length SPECIES_GAP_REGION_LENGTH]

MAF过滤工具

options:
  -h, --help            show this help message and exit
  --input INPUT
  --output OUTPUT
  --min_sequences MIN_SEQUENCES
  --required_species REQUIRED_SPECIES
  --ref_species REF_SPECIES
  --min_ref_length MIN_REF_LENGTH
  --max_gap_ratio MAX_GAP_RATIO
  --gap_region_length GAP_REGION_LENGTH
  --species_gap_ratio SPECIES_GAP_RATIO
  --species_gap_region_length SPECIES_GAP_REGION_LENGTH
```
<a name=MAF_realigner.py></a>MAF_realigner.py的参数说明
```
usage: MAF_realigner.py [-h] -i INPUT -o OUTPUT [--muscle MUSCLE] [--mafft MAFFT] [-t THRESHOLD] [-v]

MAF文件重比对工具 - 保持序列顺序版本

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        输入MAF文件
  -o OUTPUT, --output OUTPUT
                        输出MAF文件
  --muscle MUSCLE       MUSCLE5可执行文件路径
  --mafft MAFFT         MAFFT可执行文件路径
  -t THRESHOLD, --threshold THRESHOLD
                        序列长度阈值，大于此值使用MAFFT，否则使用MUSCLE5 (默认: 1000)
  -v, --verbose         详细输出
```
<a name=sv_Calling.py></a>sv_Calling.py的参数说明
```
usage: new_bao.py [-h] --maf MAF --target TARGET --reference REFERENCE [--non_target NON_TARGET] --output_prefix OUTPUT_PREFIX [--max_target_gap MAX_TARGET_GAP] [--min_ref_gap MIN_REF_GAP]
                  [--min_target_gap MIN_TARGET_GAP] [--max_ref_gap MAX_REF_GAP] [--max_gap MAX_GAP] [--inv_more_block_max INV_MORE_BLOCK_MAX] [--min_intra_size MIN_INTRA_SIZE]
                  [--max_gap_ratio MAX_GAP_RATIO]

检测目标物种的结构变异(SV)

options:
  -h, --help            show this help message and exit
  --maf MAF             输入MAF文件
  --target TARGET       目标物种名称
  --reference REFERENCE
                        参考物种名称
  --non_target NON_TARGET
                        其他非目标物种（逗号分隔）
  --output_prefix OUTPUT_PREFIX
                        输出文件前缀
  --max_target_gap MAX_TARGET_GAP
                        目标物种最大连续间隔，在缺失检测中使用
  --min_ref_gap MIN_REF_GAP
                        参考物种最小SV间隔
  --min_target_gap MIN_TARGET_GAP
                        目标物种最小SV间隔
  --max_ref_gap MAX_REF_GAP
                        参考物种最大连续间隔，在插入检测中使用
  --max_gap MAX_GAP     连续区域最大间隔，在易位和倒位检测中使用
  --inv_more_block_max INV_MORE_BLOCK_MAX
                        多块倒位检测中，所有物种允许的最大间隔阈值
  --min_intra_size MIN_INTRA_SIZE
                        块内SV的最小尺寸
  --max_gap_ratio MAX_GAP_RATIO
                        缺失检测中非目标物种允许的最大gap比例
```

<a name=conservation_calculator.py></a>conservation_calculator.py的参数说明
```
usage: conservation_calculator.py [-h] -m MAF -b BED -o OUTPUT [-r REF_SPECIES] [-v] [--version]

计算缺失变异两端序列的保守性指标

options:
  -h, --help            show this help message and exit
  -m MAF, --maf MAF     输入MAF文件（多序列比对文件）
  -b BED, --bed BED     输入BED文件（变异文件）
  -o OUTPUT, --output OUTPUT
                        输出BED文件（添加保守性指标后的结果）
  -r REF_SPECIES, --ref-species REF_SPECIES
                        参考物种名称（默认: mu_mu）
  -v, --verbose         显示详细处理信息
  --version             显示版本信息

使用示例:
  conservation_calculator.py -m alignment.maf -b variants.bed -o result.bed
  conservation_calculator.py -m alignment.maf -b variants.bed -o result.bed -r human
  conservation_calculator.py -m alignment.maf -b variants.bed -o result.bed -v
  
输入要求:
  - MAF文件: 多序列比对文件，包含a score=和s行
  - BED文件: 包含变异信息的BED格式文件，需要有特定的列结构
```

<a name=bed格式的SV文件></a>bed格式的SV文件说明
```
##我们生成共计22列的变异信息，详细记录所检测到的目标物种特异结构变异在本身基因组以及对应参考基因组上坐标信息，为用户提供多样化的候选SV筛选标准。
Ref_Chrom: 参考基因组染色体编号
Ref_Start: 参考基因组起始位置
Ref_End: 参考基因组终止位置
Ref_Strand: 参考基因组链方向信息
Target_Chrom: 目标基因组染色体编号
Target_Start: 目标基因组起始位置
Target_End: 目标基因组终止位置
Target_Strand: 目标基因组链方向信息
Original_tra_Chrom: 易位片段染色体编号，只有易位变异才存在
Original_tra_Start: 易位片段基因组起始位置，只有易位变异才存在
Original_tra_End: 易位片段基因组终止位置，只有易位变异才存在
Original_tra_Strand: 易位片段链方向信息，只有易位变异才存在
SV_Type: 基因组结构变异的类型
SV_Size: 基因组结构变异的大小(bp)
SV_Source: 变异所在比对块中的情况，比如在比对块内还是比对块间
Block_Indexes: 变异所涉及的比对块索引
Additional_Info: 额外信息，同SV_Source
SV_Upstream_Sequence_Consistency: 变异位置上游序列一致性
SV_Downstream_Sequence_Consistency: 变异位置下游序列一致性
SV_Total_Sequence_Consistency: 变异位置上下游总的序列一致性
SV_Upstream_Length: 变异位置上游序列一致性计算范围
SV_Downstream_Length: 变异位置下游序列一致性计算范围
```

## Citation:

