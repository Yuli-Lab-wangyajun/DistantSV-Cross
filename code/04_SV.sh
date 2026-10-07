#!/bin/bash

# SV Analysis Script Generator
# Usage: ./04_SV.sh -s species.txt -r reference -d working_directory -t Template3_SV_script.sh -target target_species -non_target "non_target_list"

usage() {
    echo "Usage: $0 -s species.txt -r reference_genome -d working_directory -t template.sh -target target_species -non_target \"non_target_list\""
    echo "Parameter description:"
    echo "  -s species.txt        : File containing genome abbreviations and their paths"
    echo "  -r reference_genome   : Abbreviation of the reference genome"
    echo "  -d working_directory  : Working directory path"
    echo "  -t template.sh        : SV analysis template script"
    echo "  -target target_species: Abbreviation of the target species"
    echo "  -non_target \"list\"    : Non-target species list, comma-separated"
    exit 1
}

# Manually parse parameters
SPECIES_FILE=""
REF_GENOME=""
WORK_DIR=""
TEMPLATE_FILE=""
TARGET_SPECIES=""
NON_TARGET_SPECIFIED=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -s)
            SPECIES_FILE="$2"
            shift 2
            ;;
        -r)
            REF_GENOME="$2"
            shift 2
            ;;
        -d)
            WORK_DIR="$2"
            shift 2
            ;;
        -t)
            TEMPLATE_FILE="$2"
            shift 2
            ;;
        -target)
            TARGET_SPECIES="$2"
            shift 2
            ;;
        -non_target)
            NON_TARGET_SPECIFIED="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            usage
            ;;
    esac
done

# Check necessary parameters
if [[ -z "$SPECIES_FILE" || -z "$REF_GENOME" || -z "$WORK_DIR" || -z "$TEMPLATE_FILE" || -z "$TARGET_SPECIES" || -z "$NON_TARGET_SPECIFIED" ]]; then
    echo "Error: Missing required parameter!"
    usage
fi

# Check if files exist
if [[ ! -f "$SPECIES_FILE" ]]; then
    echo "Error: species.txt file does not exist: $SPECIES_FILE"
    exit 1
fi

if [[ ! -f "$TEMPLATE_FILE" ]]; then
    echo "Error: Template file does not exist: $TEMPLATE_FILE"
    exit 1
fi

# Enter the working directory
cd "$WORK_DIR" || exit 1

# Directly use the user-specified non_target list
NON_TARGET_STRING="$NON_TARGET_SPECIFIED"
echo "Using user-specified non_target list: $NON_TARGET_STRING"

# Retrieve all chromosomes
echo "Reading chromosome information..."
CHR_MULTIZ_DIRS=()
for dir in *_multiz_out; do
    if [[ -d "$dir" ]]; then
        chr_name=$(basename "$dir" | sed 's/_multiz_out$//')
        CHR_MULTIZ_DIRS+=("$chr_name")
        echo "  Found chromosome: $chr_name"
    fi
done

# If none found above, try to infer from axt directory
if [[ ${#CHR_MULTIZ_DIRS[@]} -eq 0 ]]; then
    echo "Inferring chromosome information from axt directory..."
    for dir in axt/*_axt; do
        if [[ -d "$dir" ]]; then
            chr_name=$(basename "$dir" | sed 's/_axt$//')
            CHR_MULTIZ_DIRS+=("$chr_name")
            echo "  Inferred chromosome: $chr_name"
        fi
    done
fi

# Generate SV analysis scripts
echo "Generating SV analysis scripts..."
for chr_name in "${CHR_MULTIZ_DIRS[@]}"; do
    script_name="shell_sh/${chr_name}_SV.sh"
    
    # Read the template file
    template_content=$(cat "$TEMPLATE_FILE")
    
    # Replace job name
    new_job_name="${chr_name}_SV"
    template_content=$(echo "$template_content" | sed "s|JOB_NAME_PLACEHOLDER|$new_job_name|")
    
    # Replace chromosome name
    template_content=$(echo "$template_content" | sed "s|CHROMOSOME_PLACEHOLDER|$chr_name|g")
    
    # Replace reference genome name
    template_content=$(echo "$template_content" | sed "s|REF_GENOME_PLACEHOLDER|$REF_GENOME|g")
    
    # Replace working directory path
    template_content=$(echo "$template_content" | sed "s|WORK_DIR_PLACEHOLDER|$WORK_DIR|g")
    
    # Replace non_target list
    template_content=$(echo "$template_content" | sed "s|NON_TARGET_PLACEHOLDER|$NON_TARGET_STRING|g")
    
    # Replace target species
    template_content=$(echo "$template_content" | sed "s|TARGET_SPECIES_PLACEHOLDER|$TARGET_SPECIES|g")
    
    # Write to script file
    echo "$template_content" > "$script_name"
    chmod +x "$script_name"
    
    echo "  Generated script: $script_name"
done

# Generate batch submission script
echo "Generating batch submission script..."
BATCH_SV_SUBMIT="shell_sh/batch_sv_submit.sh"
cat > "$BATCH_SV_SUBMIT" << EOF
#!/bin/bash
# Batch submission script for SV analysis jobs

echo "Starting to submit SV analysis jobs..."

EOF

for chr_name in "${CHR_MULTIZ_DIRS[@]}"; do
    script_name="${chr_name}_SV.sh"
    echo "echo \"Submitting job: $script_name\"" >> "$BATCH_SV_SUBMIT"
    echo "sbatch $script_name" >> "$BATCH_SV_SUBMIT"
    echo "sleep 2" >> "$BATCH_SV_SUBMIT"  # Avoid submitting too many jobs simultaneously
done

cat >> "$BATCH_SV_SUBMIT" << EOF

echo "All SV analysis jobs have been submitted!"
EOF

chmod +x "$BATCH_SV_SUBMIT"

echo ""
echo "SV analysis script generation completed!"
echo "Number of generated scripts: ${#CHR_MULTIZ_DIRS[@]}"
echo "Reference genome: $REF_GENOME"
echo "Target species: $TARGET_SPECIES"
echo "Non-target list: $NON_TARGET_STRING"
echo ""
echo "Usage:"
echo "1. Submit jobs individually: cd shell_sh && sbatch [script_name]"
echo "2. Submit jobs in batch: cd shell_sh && ./batch_sv_submit.sh"
echo ""
echo "Generated scripts are located in: shell_sh/"