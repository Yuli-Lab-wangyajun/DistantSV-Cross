#!/usr/bin/env python3
"""
Extract coordinate information for specified sequences from MAF format files (positive strand only)
"""

import argparse
import sys

def parse_maf_file(maf_file, target_species):
    """
    Parse MAF file and extract coordinate information of the target species (positive strand only)
    
    Parameters:
        maf_file: Path to the MAF file
        target_species: Name of the target species
    
    Returns:
        List containing coordinate information
    """
    coordinates = []
    
    try:
        with open(maf_file, 'r') as f:
            current_block = []
            in_block = False
            
            for line in f:
                line = line.strip()
                
                if line.startswith('a'):
                    # Start of a new alignment block (compatible with lines with parameters)
                    if in_block and current_block:
                        process_block(current_block, target_species, coordinates)
                    current_block = []
                    in_block = True
                elif line.startswith('s') and in_block:
                    # Sequence line
                    current_block.append(line)
                elif not line and in_block:
                    # Empty line, end of block
                    if current_block:
                        process_block(current_block, target_species, coordinates)
                    current_block = []
                    in_block = False
            
            # Process the last block
            if in_block and current_block:
                process_block(current_block, target_species, coordinates)
                
    except FileNotFoundError:
        print(f"Error: File {maf_file} not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Error occurred while reading file - {e}", file=sys.stderr)
        sys.exit(1)
    
    return coordinates

def process_block(block_lines, target_species, coordinates):
    """
    Process a single alignment block and extract positive strand information of the target species
    """
    for line in block_lines:
        if line.startswith('s'):
            parts = line.split()
            if len(parts) >= 7:
                seq_name = parts[1]
                
                # Check if it matches the target species (species name before the dot)
                if seq_name.split('.')[0] == target_species:
                    try:
                        start = int(parts[2])
                        length = int(parts[3])
                        strand = parts[4]
                        total_length = int(parts[5])
                        sequence_data = parts[6]
                        
                        # ---------- Key modification: Extract positive strand (+) alignments only ----------
                        if strand == '+':
                            coordinates.append({
                                'species': target_species,
                                'chromosome': '.'.join(seq_name.split('.')[1:]),  # Part after the dot as chromosome
                                'start': start,
                                'length': length,
                                'end': start + length,
                                'strand': strand,
                                'total_length': total_length,
                                'sequence': sequence_data
                            })
                    except ValueError as e:
                        print(f"Warning: Error parsing coordinates - {e}", file=sys.stderr)
                        continue

def format_output(coordinates, output_format):
    """
    Format the output results
    """
    output_lines = []
    
    if output_format == 'bed':
        # BED format output
        for coord in coordinates:
            # BED format: chrom start end name score strand
            bed_line = f"{coord['chromosome']}\t{coord['start']}\t{coord['end']}\t{coord['species']}\t0\t{coord['strand']}"
            output_lines.append(bed_line)
    
    elif output_format == 'tsv':
        # TSV table format output
        header = "Species\tChromosome\tStart\tEnd\tLength\tStrand\tTotal_Length"
        output_lines.append(header)
        for coord in coordinates:
            tsv_line = f"{coord['species']}\t{coord['chromosome']}\t{coord['start']}\t{coord['end']}\t{coord['length']}\t{coord['strand']}\t{coord['total_length']}"
            output_lines.append(tsv_line)
    
    elif output_format == 'coord':
        # Coordinate format: chromosome start length
        for coord in coordinates:
            coord_line = f"{coord['chromosome']}\t{coord['start']}\t{coord['length']}"
            output_lines.append(coord_line)
    
    elif output_format == 'full':
        # Full information format
        for coord in coordinates:
            full_line = f"{coord['species']}.{coord['chromosome']}\t{coord['start']}\t{coord['length']}\t{coord['strand']}\t{coord['total_length']}\t{coord['sequence']}"
            output_lines.append(full_line)
    
    return output_lines

def main():
    parser = argparse.ArgumentParser(
        description='Extract positive-strand alignment coordinate information for specified species from MAF files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python maf_coordinate_extractor.py -i alignment.maf -s mu_mu -o output.bed -f bed
  python maf_coordinate_extractor.py -i alignment.maf -s mu_mu -o coordinates.tsv -f tsv
  python maf_coordinate_extractor.py -i alignment.maf -s mu_mu --stdout -f coord
        '''
    )
    
    parser.add_argument('-i', '--input', required=True, help='Path to input MAF file')
    parser.add_argument('-s', '--species', required=True, help='Name of species to extract coordinates for (positive strand only)')
    parser.add_argument('-o', '--output', help='Path to output file (output to stdout if not specified)')
    parser.add_argument('-f', '--format', choices=['bed', 'tsv', 'coord', 'full'], 
                       default='coord', help='Output format: bed, tsv, coord, full (default: coord)')
    parser.add_argument('--stdout', action='store_true', help='Force output to standard output')
    
    args = parser.parse_args()
    
    # Extract coordinate information
    print(f"Extracting positive-strand coordinates for species {args.species} from {args.input}...", file=sys.stderr)
    coordinates = parse_maf_file(args.input, args.species)
    
    if not coordinates:
        print(f"Warning: No positive-strand alignments for species {args.species} found in file {args.input}", file=sys.stderr)
        sys.exit(0)
    
    print(f"Found {len(coordinates)} matching positive-strand coordinates", file=sys.stderr)
    
    # Format output
    output_lines = format_output(coordinates, args.format)
    
    # Output results
    if args.output and not args.stdout:
        try:
            with open(args.output, 'w') as f:
                for line in output_lines:
                    f.write(line + '\n')
            print(f"Results saved to {args.output}", file=sys.stderr)
        except Exception as e:
            print(f"Error: Error writing to output file - {e}", file=sys.stderr)
            sys.exit(1)
    else:
        for line in output_lines:
            print(line)

if __name__ == '__main__':
    main()