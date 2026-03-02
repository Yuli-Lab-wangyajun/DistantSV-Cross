#!/bin/bash
#SBATCH -J script_test_sv
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o /work/home/acm2g39zj0/err_out/std.out.%j
#SBATCH -e /work/home/acm2g39zj0/err_out/std.err.%j

source ~/miniconda3/bin/activate
conda activate test
##step1 Genome alignment preprocessing；第一步是生成预备文件，需要手动删掉参考物种中不用的序列文件
cd /work/home/acm2g39zj0/wangyj/SV/new_SV_method/really_data/script_test
#./01_pre.sh -s species.txt -r genome1 -d /work/home/acm2g39zj0/wangyj/SV/new_SV_method/really_data/script_test
##step2 Pairwise genome alignment；批量生成两两基因组比对脚本，将其他基因组一次比对到参考基因组的每一条染色体上，节省资源
#./02_align.sh -s species.txt -r genome1 -d /work/home/acm2g39zj0/wangyj/SV/new_SV_method/really_data/script_test -t Template1_lastz_script.sh
##step3 Genome alignment preprocessing鉴于脚本是串行运行的，用户可自行整合脚本，最终运行一个脚本即可
#./03_MSA.sh -s species.txt -r genome1 -d /work/home/acm2g39zj0/wangyj/SV/new_SV_method/really_data/script_test -t Template2_multiz_script.sh -tree "(((genome1 genome2) target) genome3)"
##step4-6 Genome alignment preprocessing同上，鉴于脚本是串行运行的，用户可自行整合脚本，最终运行一个脚本即可
./04_SV.sh -s species.txt -r genome1 -d /work/home/acm2g39zj0/wangyj/SV/new_SV_method/really_data/script_test -t Template3_SV_script.sh -target target -non_target "genome1,genome2,genome3"

