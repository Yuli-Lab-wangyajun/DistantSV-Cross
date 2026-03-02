#!/bin/bash

# Batch generation of LASTZ comparison scripts
# usage: ./generate_align_scripts.sh -s species.txt -r reference_genome -d working_directory -t lastz_template.sh
usage() {
    echo "usage: $0 -s species.txt -r reference_genome -d working_directory -t template.sh"
    echo "parameter description:"
    echo "  -s species.txt        : Files of genome abbreviations and their paths"
    echo "  -r reference_genome   : Reference"
    echo "  -d working_directory  : Working directory path"
    echo "  -t template.sh        : LASTZ template script"
    exit 1
}

# Parameter Analysis
while getopts "s:r:d:t:" opt; do
    case $opt in
        s) SPECIES_FILE="$OPTARG" ;;
        r) REF_GENOME="$OPTARG" ;;
        d) WORK_DIR="$OPTARG" ;;
        t) TEMPLATE_FILE="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check necessary parameters
if [[ -z "$SPECIES_FILE" || -z "$REF_GENOME" || -z "$WORK_DIR" || -z "$TEMPLATE_FILE" ]]; then
    echo "Error: Missing required parameter!"
    usage
fi

# Check whether the species.txt file exists
if [[ ! -f "$SPECIES_FILE" ]]; then
    echo "Error: species.txt file does not exist: $SPECIES_FILE"
    exit 1
fi

if [[ ! -f "$TEMPLATE_FILE" ]]; then
    echo "Error: template file does not exist: $TEMPLATE_FILE"
    exit 1
fi

# Enter the working directory
cd "$WORK_DIR" || exit 1

# Create the necessary directory structure
mkdir -p axt shell_sh

# Obtain all genomes other than the reference genome
echo "Read non-reference genome..."
declare -A GENOME_PATHS
NON_REF_GENOMES=()

# Reset line counter
line_num=0
while IFS=$'\t' read -r species_name genome_path || [[ -n "$species_name" ]]; do
    line_num=$((line_num+1))
    # Skip empty lines and comment lines
    [[ -z "$species_name" || "$species_name" =~ ^# ]] && continue
    # Remove possible spaces
    species_name=$(echo "$species_name" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    genome_path=$(echo "$genome_path" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    # Check whether the two fields have been successfully read
    if [[ -z "$genome_path" ]]; then
        echo "Warning: Skipping lines with formatting errors $line_num: '$species_name'"
        continue
    fi
    
    # Store all genome paths
    GENOME_PATHS["$species_name"]="$genome_path"
    
    # If it is a non-reference genome, add it to the list
    if [[ "$species_name" != "$REF_GENOME" ]]; then
        NON_REF_GENOMES+=("$species_name")
        echo "  non-reference genome: $species_name -> $genome_path"
    else
        echo "  reference genome: $species_name -> $genome_path"
    fi
done < "$SPECIES_FILE"

echo "find ${#NON_REF_GENOMES[@]} non-reference genome"

# Get all chromosome files
echo "read chromosome files..."
CHR_FILES=()
for chr_file in ref_chr/*; do
    if [[ -f "$chr_file" ]]; then
        chr_name=$(basename "$chr_file")
        CHR_FILES+=("$chr_name")
        echo "  chromosome: $chr_name"
    fi
done

# Create an axt subdirectory for each chromosome
for chr_name in "${CHR_FILES[@]}"; do
    axt_subdir="axt/${chr_name}_axt"
    mkdir -p "$axt_subdir"
    echo "  create: $axt_subdir"
done

# Generate comparison script
for chr_name in "${CHR_FILES[@]}"; do
    for genome in "${NON_REF_GENOMES[@]}"; do
        # Script file name
        script_name="shell_sh/${chr_name}_${genome}_lastz.sh"
        
        # Read template file
        template_content=$(cat "$TEMPLATE_FILE")
        
        # Replace assignment name
        new_job_name="${chr_name}_${genome}_lastz"
        template_content=$(echo "$template_content" | sed "s|JOB_NAME_PLACEHOLDER|$new_job_name|")
        
        # Replace chromosome name
        template_content=$(echo "$template_content" | sed "s|CHROMOSOME_PLACEHOLDER|$chr_name|g")
        
        # Replace genome path
        template_content=$(echo "$template_content" | sed "s|GENOME_NAME_PLACEHOLDER|${GENOME_PATHS[$genome]}|g")
        
        # Replace output file name
        template_content=$(echo "$template_content" | sed "s|GENOME_NAME_OUTPUT|$genome|g")
        
        # Replace work directory path
        template_content=$(echo "$template_content" | sed "s|WORK_DIR_PLACEHOLDER|$WORK_DIR|g")
        
        # Write to script file
        echo "$template_content" > "$script_name"
        chmod +x "$script_name"
        
        echo "  Generate script: $script_name"
    done
done

# Generate batch submission script
BATCH_SUBMIT="shell_sh/batch_lastz_submit.sh"
cat > "$BATCH_SUBMIT" << EOF
#!/bin/bash
# Batch submission comparison job script

echo "Start submitting tasks..."

EOF

for chr_name in "${CHR_FILES[@]}"; do
    for genome in "${NON_REF_GENOMES[@]}"; do
        script_name="${chr_name}_${genome}_lastz.sh"
        echo "echo \"submitting tasks: $script_name\"" >> "$BATCH_SUBMIT"
        echo "sbatch $script_name" >> "$BATCH_SUBMIT"
        echo "sleep 1" >> "$BATCH_SUBMIT"  # Avoid submitting too many assignments at the same time
    done
done

cat >> "$BATCH_SUBMIT" << EOF

EOF

chmod +x "$BATCH_SUBMIT"

echo ""
echo "Alignment scripts generation completed!"