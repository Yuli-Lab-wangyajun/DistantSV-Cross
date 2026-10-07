#!/bin/bash
#SBATCH -J JOB_NAME_PLACEHOLDER
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=TASKS_PER_NODE_PLACEHOLDER
#SBATCH -o /work/home/acm2g39zj0/err_out/std.out.%j
#SBATCH -e /work/home/acm2g39zj0/err_out/std.err.%j

set -euo pipefail

conda activate your_env

WORK_DIR="WORK_DIR_PLACEHOLDER"
CHROMOSOME="CHROMOSOME_PLACEHOLDER"
GENOME_PATH="GENOME_NAME_PLACEHOLDER"
JOB_NAME="JOB_NAME_PLACEHOLDER"
OUTPUT_NAME="GENOME_NAME_OUTPUT"

cd "$WORK_DIR/axt/${CHROMOSOME}_axt"

##Calculate the length of each chromosome of the reference species
REF_LEN=$(
    awk '
        !/^>/ {
            gsub(/[[:space:]\r]/, "", $0)
            total += length($0)
        }
        END {
            print total + 0
        }
    ' "$WORK_DIR/ref_chr/$CHROMOSOME"
)
if [[ ! "$REF_LEN" =~ ^[0-9]+$ ]] || (( REF_LEN <= 0 )); then
    echo "Error: Cannot calculate chromosome length."
    exit 1
fi

##Different segments of chromosomes of different lengths run in parallel
if (( REF_LEN >= 150000000 )); then
    SEGMENT_COUNT=20
elif (( REF_LEN >= 100000000 )); then
    SEGMENT_COUNT=10
elif (( REF_LEN >= 30000000 )); then
    SEGMENT_COUNT=5
else
    SEGMENT_COUNT=1
fi

TMP_DIR="tmp_chunk_dir_${JOB_NAME}"
mkdir -p "$TMP_DIR"
rm -f "$TMP_DIR"/c*.axt

##Lastz Running Module
run_lastz() {
    local chunk_id="$1"
    local start="$2"
    local end="$3"
    echo "Start chunk ${chunk_id}: ${start}..${end}"
    lastz_32 \
        "$WORK_DIR/ref_chr/$CHROMOSOME[${start}..${end}]" \
        "$GENOME_PATH" \
        --ambiguous=iupac \
        --inner=2000 \
        --step=1 \
        --ydrop=3000 \
        --gappedthresh=6000 \
        --hspthresh=2200 \
        --allocate:traceback=2G \
        --format=axt \
        --out="$TMP_DIR/c${chunk_id}.axt"
}
PIDS=()

##Uniformly divide the chromosomes and execute in parallel
for ((i=1; i<=SEGMENT_COUNT; i++)); do
    START=$(( ((i - 1) * REF_LEN) / SEGMENT_COUNT + 1 ))
    END=$(( (i * REF_LEN) / SEGMENT_COUNT ))
    if (( i == SEGMENT_COUNT )); then
        END="$REF_LEN"
    fi
    run_lastz "$i" "$START" "$END" &
    PIDS+=("$!")
done
FAILED=0
for pid in "${PIDS[@]}"; do
    if ! wait "$pid"; then
        FAILED=1
    fi
done
if (( FAILED != 0 )); then
    echo "Error: At least one LASTZ chunk failed."
    exit 1
fi

##Merge the axt files
python3 - "$TMP_DIR" "$OUTPUT_NAME" "$SEGMENT_COUNT" <<'PY'
import os
import re
import sys

tmp_dir = sys.argv[1]
output_name = sys.argv[2]
segment_count = int(sys.argv[3])

subrange_pattern = re.compile(r"\[\d+\.\.\d+\]")
block_id = 0

with open(f"{output_name}.axt", "w") as out:
    for i in range(1, segment_count + 1):
        file_path = os.path.join(tmp_dir, f"c{i}.axt")
        if not os.path.exists(file_path):
            raise FileNotFoundError(
                f"Missing LASTZ chunk output: {file_path}"
            )
        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) == 9 and parts[0].isdigit():
                    parts[0] = str(block_id)
                    block_id += 1
                    parts[1] = subrange_pattern.sub("", parts[1])
                    out.write(" ".join(parts) + "\n")
                else:
                    out.write(line)
PY

##Delete temporary files
rm -rf "$TMP_DIR"