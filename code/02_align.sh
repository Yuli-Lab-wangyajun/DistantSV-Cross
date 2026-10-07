#!/bin/bash

# Batch generation of LASTZ comparison scripts
set -euo pipefail

usage() {
    echo "usage: $0 -s species.txt -r reference_genome -d working_directory -t template.sh"
    echo ""
    echo "parameter description:"
    echo "  -s species.txt        : Files of genome abbreviations and their paths"
    echo "  -r reference_genome   : Reference genome name"
    echo "  -d working_directory  : Working directory path"
    echo "  -t template.sh        : LASTZ template script"
    exit 1
}

SPECIES_FILE=""
REF_GENOME=""
WORK_DIR=""
TEMPLATE_FILE=""

while getopts "s:r:d:t:" opt; do
    case $opt in
        s) SPECIES_FILE="$OPTARG" ;;
        r) REF_GENOME="$OPTARG" ;;
        d) WORK_DIR="$OPTARG" ;;
        t) TEMPLATE_FILE="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$SPECIES_FILE" ||
      -z "$REF_GENOME" ||
      -z "$WORK_DIR" ||
      -z "$TEMPLATE_FILE" ]]; then
    echo "Error: Missing required parameter!"
    usage
fi

if [[ ! -f "$SPECIES_FILE" ]]; then
    echo "Error: species.txt file does not exist: $SPECIES_FILE"
    exit 1
fi

if [[ ! -f "$TEMPLATE_FILE" ]]; then
    echo "Error: template file does not exist: $TEMPLATE_FILE"
    exit 1
fi

SPECIES_FILE="$(readlink -f "$SPECIES_FILE")"
TEMPLATE_FILE="$(readlink -f "$TEMPLATE_FILE")"
WORK_DIR="$(readlink -f "$WORK_DIR")"

cd "$WORK_DIR" || exit 1

mkdir -p axt shell_sh

escape_sed_value() {
    printf '%s' "$1" | sed 's/[\/&|\\]/\\&/g'
}

replace_placeholder() {
    local file="$1"
    local placeholder="$2"
    local value="$3"

    local escaped_value
    escaped_value="$(escape_sed_value "$value")"

    sed -i "s|${placeholder}|${escaped_value}|g" "$file"
}

get_fasta_length() {
    local fasta_file="$1"

    awk '
        !/^>/ {
            gsub(/[[:space:]\r]/, "", $0)
            total += length($0)
        }
        END {
            print total + 0
        }
    ' "$fasta_file"
}
get_segment_count() {
    local chr_length="$1"

    if (( chr_length >= 150000000 )); then
        echo 20
    elif (( chr_length >= 100000000 )); then
        echo 10
    elif (( chr_length >= 30000000 )); then
        echo 5
    else
        echo 1
    fi
}

echo "Read non-reference genome..."

declare -A GENOME_PATHS
NON_REF_GENOMES=()

line_num=0

while IFS=$'\t' read -r species_name genome_path || [[ -n "$species_name" ]]; do
    line_num=$((line_num + 1))
    species_name="${species_name//$'\r'/}"
    genome_path="${genome_path//$'\r'/}"
    species_name="$(echo "$species_name" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    genome_path="$(echo "$genome_path" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    [[ -z "$species_name" ]] && continue
    [[ "$species_name" =~ ^# ]] && continue
    if [[ -z "$genome_path" ]]; then
        echo "Warning: Skipping malformed line ${line_num}: $species_name"
        continue
    fi
    GENOME_PATHS["$species_name"]="$genome_path"

    if [[ "$species_name" != "$REF_GENOME" ]]; then
        NON_REF_GENOMES+=("$species_name")
        echo "  non-reference genome: $species_name -> $genome_path"
    else
        echo "  reference genome: $species_name -> $genome_path"
    fi
done < "$SPECIES_FILE"

echo "find ${#NON_REF_GENOMES[@]} non-reference genome(s)"

if [[ ${#NON_REF_GENOMES[@]} -eq 0 ]]; then
    echo "Error: No non-reference genome found!"
    exit 1
fi

echo "read chromosome files..."

shopt -s nullglob
CHR_FILES=(ref_chr/*)
shopt -u nullglob

CHR_FILES=(
    $(printf '%s\n' "${CHR_FILES[@]}" | sort -V)
)

if [[ ${#CHR_FILES[@]} -eq 0 ]]; then
    echo "Error: No chromosome files found in ref_chr/"
    exit 1
fi

##Create an AXT directory for each chromosome
for chr_file in "${CHR_FILES[@]}"; do
    if [[ -f "$chr_file" ]]; then
        chr_name="$(basename "$chr_file")"
        mkdir -p "axt/${chr_name}_axt"
        echo "  chromosome: $chr_name"
    fi
done
MANIFEST="shell_sh/lastz_manifest.tsv"
cat > "$MANIFEST" <<EOF
EOF
echo ""
echo "Generating LASTZ comparison scripts..."

for chr_file in "${CHR_FILES[@]}"; do
    [[ -f "$chr_file" ]] || continue
    chr_name="$(basename "$chr_file")"
    chr_length="$(get_fasta_length "$chr_file")"
    if [[ ! "$chr_length" =~ ^[0-9]+$ ]] || (( chr_length <= 0 )); then
        echo "Warning: Invalid chromosome length, skip: $chr_file"
        continue
    fi
    segment_count="$(get_segment_count "$chr_length")"
    requested_cores="$segment_count"
    echo ""
    echo "Chromosome: $chr_name"
    echo "  length : $chr_length bp"
    echo "  chunks : $segment_count"
    echo "  cores  : $requested_cores"
    for genome in "${NON_REF_GENOMES[@]}"; do
        genome_path="${GENOME_PATHS[$genome]}"
        if [[ "$genome_path" != /* ]]; then
            genome_path="$WORK_DIR/$genome_path"
        fi
        script_name="shell_sh/${chr_name}_${genome}_lastz.sh"
        job_name="${chr_name}_${genome}_lastz"
        cp "$TEMPLATE_FILE" "$script_name"
        replace_placeholder \
            "$script_name" \
            "JOB_NAME_PLACEHOLDER" \
            "$job_name"
        replace_placeholder \
            "$script_name" \
            "CHROMOSOME_PLACEHOLDER" \
            "$chr_name"
        replace_placeholder \
            "$script_name" \
            "GENOME_NAME_PLACEHOLDER" \
            "$genome_path"
        replace_placeholder \
            "$script_name" \
            "GENOME_NAME_OUTPUT" \
            "$genome"
        replace_placeholder \
            "$script_name" \
            "WORK_DIR_PLACEHOLDER" \
            "$WORK_DIR"
        replace_placeholder \
            "$script_name" \
            "TASKS_PER_NODE_PLACEHOLDER" \
            "$requested_cores"
        chmod +x "$script_name"
        absolute_script_name="$(readlink -f "$script_name")"
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$genome" \
            "$chr_name" \
            "$chr_length" \
            "$segment_count" \
            "$requested_cores" \
            "$absolute_script_name" \
            >> "$MANIFEST"
        echo "  Generate script: $script_name"
    done
done