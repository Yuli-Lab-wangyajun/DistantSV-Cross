#!/bin/bash

# Batch generation of MULTIZ alignment scripts
# usage: ./MSA.sh -s species.txt -r reference -d working_directory -t multiz_template.sh -tree "species_tree"

usage() {
    echo "usage: $0 -s species.txt -r reference_genome -d working_directory -t multiz_template.sh -tree \"species_tree\""
    echo "parameter description:"
    echo "  -s species.txt        : Files of genome abbreviations and their paths"
    echo "  -r reference_genome   : Reference"
    echo "  -d working_directory  : Working directory path"
    echo "  -t template.sh        : MULTIZ template script"
    echo "  -tree \"species_tree\"  : species tree (required)"
    exit 1
}

# Parameter Analysis
SPECIES_FILE=""
REF_GENOME=""
WORK_DIR=""
TEMPLATE_FILE=""
SPECIES_TREE=""

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
        -tree)
            SPECIES_TREE="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            usage
            ;;
    esac
done

# Check necessary parameters
if [[ -z "$SPECIES_FILE" || -z "$REF_GENOME" || -z "$WORK_DIR" || -z "$TEMPLATE_FILE" || -z "$SPECIES_TREE" ]]; then
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

# Check if necessary directories exist
if [[ ! -d "axt" || ! -d "species_nibs" || ! -d "species_sizes" || ! -d "species_2bit" ]]; then
    echo "Error: Necessary directories are missing. Please run previous scripts first to generate required files."
    exit 1
fi

# Retrieve all chromosomes
echo "Reading chromosome information..."
CHR_AXT_DIRS=()
for dir in axt/*_axt; do
    if [[ -d "$dir" ]]; then
        chr_name=$(basename "$dir" | sed 's/_axt$//')
        CHR_AXT_DIRS+=("$chr_name")
        echo "Found chromosome: $chr_name"
    fi
done

# Generate multiple sequence alignment scripts
echo "Generating multiple sequence alignment scripts..."
for chr_name in "${CHR_AXT_DIRS[@]}"; do
    script_name="shell_sh/${chr_name}_multiz.sh"
    
    # Read the template file
    template_content=$(cat "$TEMPLATE_FILE")
    
    # Replace job name
    new_job_name="${chr_name}_multiz"
    template_content=$(echo "$template_content" | sed "s|JOB_NAME_PLACEHOLDER|$new_job_name|")
    
    # Replace output and error log paths
    template_content=$(echo "$template_content" | sed "s|multiz_CHROMOSOME_PLACEHOLDER.out|multiz_${chr_name}.out|g")
    template_content=$(echo "$template_content" | sed "s|multiz_CHROMOSOME_PLACEHOLDER.err|multiz_${chr_name}.err|g")
    
    # Replace chromosome name
    template_content=$(echo "$template_content" | sed "s|CHROMOSOME_PLACEHOLDER|$chr_name|g")
    
    # Replace working directory path
    template_content=$(echo "$template_content" | sed "s|WORK_DIR_PLACEHOLDER|$WORK_DIR|g")
    
    # Replace reference genome name
    template_content=$(echo "$template_content" | sed "s|REF_GENOME_PLACEHOLDER|$REF_GENOME|g")
    
    # Replace species tree
    template_content=$(echo "$template_content" | sed "s|SPECIES_TREE_PLACEHOLDER|$SPECIES_TREE|")
    
    # Write to script file
    echo "$template_content" > "$script_name"
    chmod +x "$script_name"
    
    echo "Generated script: $script_name"
done

# Generate batch submission script
echo "Generating batch submission script..."
BATCH_MULTIZ_SUBMIT="shell_sh/batch_multiz_submit.sh"
cat > "$BATCH_MULTIZ_SUBMIT" << EOF
#!/bin/bash
# Batch submission script for multiple sequence alignment jobs

echo "Starting to submit multiple sequence alignment jobs..."

EOF

for chr_name in "${CHR_AXT_DIRS[@]}"; do
    script_name="${chr_name}_multiz.sh"
    echo "echo \"Submitting job: $script_name\"" >> "$BATCH_MULTIZ_SUBMIT"
    echo "sbatch $script_name" >> "$BATCH_MULTIZ_SUBMIT"
    echo "sleep 2" >> "$BATCH_MULTIZ_SUBMIT"  # Avoid submitting too many jobs simultaneously
done

cat >> "$BATCH_MULTIZ_SUBMIT" << EOF

echo "All multiple sequence alignment jobs have been submitted!"
EOF

chmod +x "$BATCH_MULTIZ_SUBMIT"

echo ""
echo "Multiple sequence alignment script generation completed!"
echo "Number of generated scripts: ${#CHR_AXT_DIRS[@]}"
echo "Species tree used: $SPECIES_TREE"
echo ""
echo "Usage:"
echo "1. Submit jobs individually: cd shell_sh && sbatch [script_name]"
echo "2. Submit jobs in batch: cd shell_sh && ./batch_multiz_submit.sh"
echo ""
echo "Generated scripts are located in: shell_sh/"
echo "Multiple sequence alignment results will be output to: [chromosome_name_multiz_out]/"