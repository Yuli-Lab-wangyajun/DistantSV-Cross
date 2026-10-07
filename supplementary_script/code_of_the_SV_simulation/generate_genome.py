import argparse
import random
import os

def generate_chromosome(output_file, species_name, chr_num, size):
    """Generate a single chromosome FASTA sequence and write to file"""
    output_file.write(f'>{species_name}_chr{chr_num}\n')  # Use species name extracted from filename
    remaining = size
    while remaining > 0:
        line_length = min(80, remaining)
        line = ''.join(random.choices(['A','T','C','G'], k=line_length))
        output_file.write(line + '\n')
        remaining -= line_length

def extract_species_name(output_path):
    """Extract species name from output file path"""
    # Get filename (without path)
    filename = os.path.basename(output_path)
    # Remove extension (e.g., .fa, .fasta)
    species_name = os.path.splitext(filename)[0]
    return species_name

def main():
    parser = argparse.ArgumentParser(description='Generate a FASTA file containing five chromosome sequences of specified size')
    parser.add_argument('-o', '--output', required=True, help='Output file name')
    parser.add_argument('-s', '--size', type=int, required=True, 
                       help='Base count for each chromosome (multiples of 80 recommended)')
    
    args = parser.parse_args()
    
    # Extract species name from output filename
    species_name = extract_species_name(args.output)
    
    with open(args.output, 'w') as f:
        for i in range(1, 6):  # Generate chr1 to chr5
            generate_chromosome(f, species_name, i, args.size)
    
    print(f"Successfully generated genome file: {args.output}")
    print(f"Species name: {species_name}")
    print(f"Each chromosome size: {args.size:,} bp")
    print(f"Chromosome count: 5")

if __name__ == '__main__':
    main()