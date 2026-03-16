import random
import argparse
import sys
import os
from datetime import datetime
import time

def generate_homologous_sequences(length, divergence):
    """
    Generate a pair of homologous DNA sequences with a specified degree of divergence
    """
    bases = ['A', 'T', 'C', 'G']
    
    # Generate a random ancestral sequence
    ancestor = ''.join(random.choices(bases, k=length))
    
    num_mutations = int(divergence * length)
    
    # Ensure that mutation positions are completely random
    mutation_positions = random.sample(range(length), num_mutations)
    
    seq1 = list(ancestor)
    seq2 = list(ancestor)
    
    for pos in mutation_positions:
        original_base = ancestor[pos]
        new_bases = [b for b in bases if b != original_base]
        new_base = random.choice(new_bases)
        
        # randomly assigned to one of the subsequences
        if random.random() < 0.5:
            seq1[pos] = new_base
        else:
            seq2[pos] = new_base
    
    return ''.join(seq1), ''.join(seq2)

def calculate_divergence(seq1, seq2):
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    return mismatches / len(seq1)

def format_fasta(sequence, header, width=80):
    formatted = [f">{header}"]
    for i in range(0, len(sequence), width):
        formatted.append(sequence[i:i+width])
    return "\n".join(formatted)

def main():
    parser = argparse.ArgumentParser(
        description='Generate completely independent homologous DNA sequence pairs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-l', '--length', type=int, default=1000,
                        help='Sequence length (bp)')
    parser.add_argument('-d', '--divergence', type=float, default=0.2,
                        help='Target Deviation (0-1)')
    parser.add_argument('-n', '--num_pairs', type=int, default=1,
                        help='Number of sequence pairs to generate')
    parser.add_argument('-o', '--output', default='homologous_sequences.fa',
                        help='Output filename format (using {i} as the sequence number)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show details')
    parser.add_argument('--id_prefix', default='seq',
                        help='Sequence ID prefix (auto-increment)')
    
    args = parser.parse_args()
    
    # Validate parameters
    if not (0 <= args.divergence <= 1):
        print("Error: The similarity score must be between 0 and 1", file=sys.stderr)
        sys.exit(1)
    
    if args.length <= 0:
        print("Error: The sequence length must be a positive integer", file=sys.stderr)
        sys.exit(1)
    
    if args.num_pairs <= 0:
        print("Error: The number of sequence pairs must be a positive integer.", file=sys.stderr)
        sys.exit(1)
    
    # Prepare the output directory
    output_dir = os.path.dirname(args.output) or '.'
    os.makedirs(output_dir, exist_ok=True)
    
    # Calculate the number of digits required for the filename format
    digits = len(str(args.num_pairs))
    
    # Generate completely independent sequence pairs
    start_time = time.time()
    for i in range(1, args.num_pairs + 1):
        # Use the system time to increase randomness
        random.seed(None)
        
        # Generate a sequence
        seq1, seq2 = generate_homologous_sequences(
            length=args.length,
            divergence=args.divergence
        )
        
        # Calculate the actual variance
        actual_divergence = calculate_divergence(seq1, seq2)
        
        # sequence ID
        id1 = f"{args.id_prefix}A_{i:0{digits}d}"
        id2 = f"{args.id_prefix}B_{i:0{digits}d}"
        
        if "{i}" in args.output:
            filename = args.output.format(i=i)
        else:
            base, ext = os.path.splitext(args.output)
            filename = f"{base}_{i:0{digits}d}{ext}"
        
        fasta_content = [
            format_fasta(seq1, id1),
            format_fasta(seq2, id2)
        ]
        
        with open(filename, 'w') as f:
            f.write("\n".join(fasta_content))
#            f.write(f"\n# Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
#            f.write(f"\n# Parameters: length={args.length}, divergence={args.divergence:.4f}")
#            f.write(f"\n# Actual divergence: {actual_divergence:.6f}")
        
        if args.verbose and i % 100 == 0:
            elapsed = time.time() - start_time
            print(f"Generated {i}/{args.num_pairs} sequences, Time required: {elapsed:.2f}seconds")
    
    total_time = time.time() - start_time
    print(f"Successfully generated {args.num_pairs} pairs of completely independent sequences!")
    print(f"Output directory: {os.path.abspath(output_dir)}")
    print(f"Time required: {total_time:.2f} seconds")
    print(f"The average per pair: {total_time/args.num_pairs*1000:.2f} milliseconds")

if __name__ == "__main__":
    main()