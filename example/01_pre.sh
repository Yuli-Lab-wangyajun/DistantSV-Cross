#!/bin/bash

# usage: ./01_pre.sh -s species.txt -r reference_genome -d working_directory
usage() {
    echo "usage: $0 -s species.txt -r reference_genome -d working_directory"
    echo "parameter description:"
    echo "  -s species.txt      : Files of genome abbreviations and their paths"
    echo "  -r reference_genome   : Reference"
    echo "  -d working_directory  : Working directory path"
    exit 1
}

# Parameter Analysis
while getopts "s:r:d:" opt; do
    case $opt in
        s) SPECIES_FILE="$OPTARG" ;;
        r) REF_GENOME="$OPTARG" ;;
        d) WORK_DIR="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check necessary parameters
if [[ -z "$SPECIES_FILE" || -z "$REF_GENOME" || -z "$WORK_DIR" ]]; then
    echo "Error: Missing required parameter!"
    usage
fi

# Check whether the species.txt file exists
if [[ ! -f "$SPECIES_FILE" ]]; then
    echo "Error: species.txt file does not exist: $SPECIES_FILE"
    exit 1
fi

# Create a working directory (if it doesn't exist)
mkdir -p "$WORK_DIR"
cd "$WORK_DIR" || exit 1

# Create the necessary directory structure
mkdir -p species_2bit species_nibs species_sizes ref_chr

# Check whether the UCSC sub-tools exist
check_command() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: Command '$1' not found. Please make sure the UCSC toolkit is installed"
        exit 1
    fi
}
check_command "faSize"
check_command "faToTwoBit"
check_command "faSplit"
check_command "faToNib"

# Parse the species.txt file and process each genome
process_genome() {
    local species_name=$1
    local genome_path=$2
    
    echo "process genome: $species_name"
    
    # Check if the genome file exists
    if [[ ! -f "$genome_path" ]]; then
        echo "Warning: Genome file does not exist: $genome_path, skipping"
        return 1
    fi
    
    # 1. Generate the genome sizes file
    faSize "$genome_path" -detailed > "species_sizes/${species_name}.sizes"
    
    # 2. Generate the genome 2bit file
    faToTwoBit "$genome_path" "species_2bit/${species_name}.2bit"
    
    # 3. Generate the genome nib file
    local nib_dir="species_nibs/${species_name}.nib"
    mkdir -p "$nib_dir"
    
    # Split genome file
    faSplit byName "$genome_path" "./$nib_dir/"
    
    # Go into the nib directory and process each file
    cd "$nib_dir" || return 1
    for file in *.fa; do
        if [[ -f "$file" ]]; then
            # Remove the .fa extension
            local base_name="${file%.fa}"
            # Convert to nib format
            faToNib -softMask "$file" "${base_name}.nib"
            # Do not delete the original fasta file anymore; keep both fa and nib files.
            echo "    Generate: ${base_name}.nib (Reserved ${file})"
        fi
    done
    
    cd - > /dev/null || return 1
    echo "  Complete the processing of $species_name"
}

# Main processing procedure
echo "Starting to process the genomes in species.txt..."
while IFS=$'\t' read -r species_name genome_path || [[ -n "$species_name" ]]; do
    # Skip empty lines and comment lines
    [[ -z "$species_name" || "$species_name" =~ ^# ]] && continue
    # Remove possible spaces
    species_name=$(echo "$species_name" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    genome_path=$(echo "$genome_path" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    # Check whether the two fields have been successfully read
    if [[ -z "$genome_path" ]]; then
        echo "Warning: Skipping lines with formatting errors: $species_name"
        continue
    fi
    
    # Processing the genome
    process_genome "$species_name" "$genome_path"
    
    # If it is a reference genome, additionally process chromosome splitting
    if [[ "$species_name" == "$REF_GENOME" ]]; then
        # Clear the ref_chr directory (if any files exist)
        rm -rf ref_chr/*
        cd ref_chr || exit 1
        faSplit byName "$genome_path" ./
        
        # File name：Remove sequence descriptions (content after the space) and the .fa extension
        for file in *; do
            if [[ -f "$file" ]]; then
                # Get the part before the first space as the new file name
                new_name=$(echo "$file" | cut -d' ' -f1)
                # Remove the .fa suffix (if any)
                new_name="${new_name%.fa}"
                if [[ "$file" != "$new_name" ]]; then
                    mv "$file" "$new_name"
                    echo "    Rename: $file -> $new_name"
                fi
            fi
        done
        
        cd - > /dev/null || exit 1
    fi
done < "$SPECIES_FILE"

echo "All processing completed!"
echo "Generated file structure:"
echo "  species_2bit/    - Genome files in 2bit format"
echo "  species_nibs/    - Genome files in nib and fa formats (split by sequence)"
echo "  species_sizes/   - Genome size information files"
echo "  ref_chr/         - Reference genome chromosome sequence files (without file extensions)"
