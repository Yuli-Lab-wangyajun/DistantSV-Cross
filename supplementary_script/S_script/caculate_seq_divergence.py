#!/usr/bin/env python3
"""
Calculate divergence between two sequences in FASTA file
"""

import argparse
import sys
from typing import Tuple

def read_fasta_sequences(fasta_file: str) -> Tuple[str, str]:
    """
    Read two sequences from FASTA file
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        Tuple containing two sequences (seq1, seq2)
    """
    sequences = []
    current_seq = []
    current_header = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Save previous sequence
                if current_header is not None and current_seq:
                    sequences.append((''.join(current_seq), current_header))
                    current_seq = []
                
                current_header = line[1:].split()[0]  # Take first word as header
            else:
                current_seq.append(line.upper())
        
        # Save last sequence
        if current_header is not None and current_seq:
            sequences.append((''.join(current_seq), current_header))
    
    if len(sequences) < 2:
        raise ValueError(f"FASTA file must contain at least 2 sequences, but found only {len(sequences)}")
    
    return sequences[0][0], sequences[1][0], sequences[0][1], sequences[1][1]

def calculate_divergence(seq1: str, seq2: str) -> float:
    """
    Calculate divergence between two sequences
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        
    Returns:
        Divergence (float between 0-1)
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequence lengths are unequal: seq1={len(seq1)}, seq2={len(seq2)}")
    
    mismatches = 0
    total_length = len(seq1)
    
    # Compare sequences base by base
    for base1, base2 in zip(seq1, seq2):
        if base1 != base2:
            mismatches += 1
    
    return mismatches / total_length

def validate_sequences(seq1: str, seq2: str) -> None:
    """
    Validate sequence validity
    """
    valid_bases = {'A', 'T', 'C', 'G', 'N'}
    
    for i, (base1, base2) in enumerate(zip(seq1, seq2)):
        if base1 not in valid_bases:
            print(f"Warning: Sequence 1 has non-standard base '{base1}' at position {i}")
        if base2 not in valid_bases:
            print(f"Warning: Sequence 2 has non-standard base '{base2}' at position {i}")

def main():
    parser = argparse.ArgumentParser(
        description='Calculate divergence between two sequences in FASTA file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  python sequence_divergence.py -i sequences.fa
  python sequence_divergence.py -i sequences.fa -v
  python sequence_divergence.py -i sequences.fa -o result.txt
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Path to input FASTA file')
    parser.add_argument('-o', '--output',
                        help='Path to output result file (optional)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show detailed information')
    
    args = parser.parse_args()
    
    try:
        # Read sequences
        if args.verbose:
            print(f"Reading FASTA file: {args.input}")
        
        seq1, seq2, header1, header2 = read_fasta_sequences(args.input)
        
        if args.verbose:
            print(f"Sequence 1: {header1}, Length: {len(seq1)}")
            print(f"Sequence 2: {header2}, Length: {len(seq2)}")
        
        # Validate sequences
        validate_sequences(seq1, seq2)
        
        # Calculate divergence
        divergence = calculate_divergence(seq1, seq2)
        
        # Calculate additional statistics
        mismatches = int(divergence * len(seq1))
        similarity = 1 - divergence
        
        # Output results
        result = f"""Sequence Divergence Analysis Result:
Input file: {args.input}
Sequence 1: {header1} (Length: {len(seq1)})
Sequence 2: {header2} (Length: {len(seq2)})
----------------------------------------
Divergence: {divergence:.6f} ({divergence * 100:.4f}%)
Similarity: {similarity:.6f} ({similarity * 100:.4f}%)
Mismatched positions: {mismatches} / {len(seq1)}
"""
        # Print to screen
        print(result)
        
        # Save to file (if specified)
        if args.output:
            with open(args.output, 'w') as f:
                f.write(result)
            if args.verbose:
                print(f"Results saved to: {args.output}")
                
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()