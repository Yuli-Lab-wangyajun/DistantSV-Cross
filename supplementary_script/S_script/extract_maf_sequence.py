#!/usr/bin/env python3
"""
Extract sequences from MAF format file and generate FASTA files
Usage: python maf_extractor.py input.maf output_prefix [--pairs PAIRS]
"""

import sys
import os
import re
import argparse
from collections import defaultdict

def parse_maf_file(maf_file):
    """
    Parse MAF file, return list of all alignment blocks
    Each block is a list of dicts containing sequence info
    """
    blocks = []
    current_block = []
    
    with open(maf_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('a'):
                # Start of new alignment block
                if current_block:
                    blocks.append(current_block)
                current_block = [{'header': line}]
            elif line.startswith('s'):
                # Sequence line
                parts = line.split()
                if len(parts) >= 7:
                    seq_info = {
                        'source': parts[1],
                        'start': int(parts[2]),
                        'size': int(parts[3]),
                        'strand': parts[4],
                        'srcSize': int(parts[5]),
                        'sequence': ''.join(parts[6:])
                    }
                    current_block.append(seq_info)
    
    # Add last block
    if current_block:
        blocks.append(current_block)
    
    return blocks

def remove_gap_columns(sequences):
    """
    Remove columns with gap ("-") in all sequences
    Return processed sequence list
    """
    if not sequences:
        return sequences
    
    # Convert sequences to char lists
    seq_chars = [list(seq['sequence']) for seq in sequences]
    seq_length = len(seq_chars[0])
    
    # Find column indices to keep (no "-" in all sequences)
    keep_columns = []
    for col in range(seq_length):
        if all(seq_chars[row][col] != '-' for row in range(len(seq_chars))):
            keep_columns.append(col)
    
    # Build new gap-free sequences
    new_sequences = []
    for i, seq in enumerate(sequences):
        new_seq = ''.join(seq_chars[i][col] for col in keep_columns)
        new_sequences.append({
            'source': seq['source'],
            'start': seq['start'],
            'size': seq['size'],
            'strand': seq['strand'],
            'srcSize': seq['srcSize'],
            'sequence': new_seq
        })
    
    print(f"    Sequence length after gap removal: {len(new_sequences[0]['sequence'])} (original: {len(sequences[0]['sequence'])})")
    return new_sequences

def concatenate_sequences_across_blocks(blocks, num_sequences=4):
    """
    Concatenate sequences at same position across all blocks
    Return concatenated sequence list
    """
    if not blocks:
        return []
    
    # Determine number of sequences to process (minimum value)
    available_sequences = min(len(block) - 1 for block in blocks)  # Minus header
    
    if available_sequences < num_sequences:
        print(f"Warning: Requested {num_sequences} sequences, but some blocks only have {available_sequences}")
        num_sequences = available_sequences
    
    # Initialize concatenated sequences
    concatenated_seqs = [{'source': '', 'sequence': ''} for _ in range(num_sequences)]
    
    for block_idx, block in enumerate(blocks):
        # Skip header, get sequence part
        sequences = block[1:1+num_sequences]
        
        # Remove gap columns
        processed_seqs = remove_gap_columns(sequences)
        
        # Concatenate sequences
        for i in range(num_sequences):
            if block_idx == 0:
                # First block, set source name
                concatenated_seqs[i]['source'] = processed_seqs[i]['source']
            concatenated_seqs[i]['sequence'] += processed_seqs[i]['sequence']
    
    return concatenated_seqs

def write_fasta_pair(seq1, seq2, filename):
    """
    Write two sequences to FASTA file
    """
    with open(filename, 'w') as f:
        # Write first sequence
        f.write(f">{seq1['source']}\n")
        f.write(f"{seq1['sequence']}\n")
        
        # Write second sequence
        f.write(f">{seq2['source']}\n")
        f.write(f"{seq2['sequence']}\n")

def parse_pair_argument(pairs_str, max_index):
    """
    Parse sequence pair argument
    Supported format: "1+2,1+3,1+4" or "all"
    """
    if pairs_str.lower() == "all":
        # Generate all possible sequence pairs
        pairs = []
        for i in range(max_index):
            for j in range(i+1, max_index):
                pairs.append((i, j))
        return pairs
    
    pairs = []
    for pair_str in pairs_str.split(','):
        try:
            a, b = pair_str.split('+')
            idx1 = int(a.strip()) - 1  # Convert to 0-based index
            idx2 = int(b.strip()) - 1
            
            if idx1 < 0 or idx1 >= max_index or idx2 < 0 or idx2 >= max_index:
                print(f"Warning: Sequence index {idx1+1}+{idx2+1} out of range, skipped")
                continue
                
            pairs.append((idx1, idx2))
        except ValueError:
            print(f"Warning: Cannot parse sequence pair '{pair_str}', skipped")
    
    return pairs

def main():
    parser = argparse.ArgumentParser(description='Extract sequences from MAF format file and generate FASTA files')
    parser.add_argument('input_maf', help='Input MAF file')
    parser.add_argument('output_prefix', help='Output file prefix')
    parser.add_argument('--pairs', default='1+2,1+3,1+4', 
                       help='Sequence pairs to generate, format like "1+2,1+3,2+3" or "all" for all combinations (default: 1+2,1+3,1+4)')
    parser.add_argument('--num_sequences', type=int, default=4,
                       help='Number of sequences to process (default: 4)')
    
    args = parser.parse_args()
    
    maf_file = args.input_maf
    output_prefix = args.output_prefix
    pairs_str = args.pairs
    num_sequences = args.num_sequences
    
    if not os.path.exists(maf_file):
        print(f"Error: Input file '{maf_file}' does not exist")
        sys.exit(1)
    
    # Parse MAF file
    print(f"Parsing MAF file: {maf_file}")
    blocks = parse_maf_file(maf_file)
    
    if not blocks:
        print("Error: No valid alignment blocks found")
        sys.exit(1)
    
    print(f"Found {len(blocks)} alignment blocks")
    
    # Check if each block has enough sequences
    available_sequences = min(len(block) - 1 for block in blocks)  # Minus header
    if available_sequences < num_sequences:
        print(f"Warning: Requested {num_sequences} sequences, but some blocks only have {available_sequences}")
        num_sequences = available_sequences
    
    print(f"Each block has at least {available_sequences} sequences, will process first {num_sequences}")
    
    # Concatenate sequences at same position across all blocks
    print("\nConcatenating sequences across all blocks and removing gaps...")
    concatenated_seqs = concatenate_sequences_across_blocks(blocks, num_sequences)
    
    if len(concatenated_seqs) < num_sequences:
        print(f"Error: Insufficient sequences ({len(concatenated_seqs)}) after concatenation (need {num_sequences})")
        sys.exit(1)
    
    print(f"\nSequence info after concatenation:")
    for i, seq in enumerate(concatenated_seqs):
        print(f"  Sequence {i+1}: {seq['source']} (length: {len(seq['sequence'])})")
    
    # Parse sequence pair argument
    pairs = parse_pair_argument(pairs_str, num_sequences)
    
    if not pairs:
        print("Error: No valid sequence pairs")
        sys.exit(1)
    
    # Create FASTA files
    print(f"\nGenerating FASTA files:")
    for idx1, idx2 in pairs:
        filename = f"{output_prefix}_seq{idx1+1}+{idx2+1}.fasta"
        write_fasta_pair(concatenated_seqs[idx1], concatenated_seqs[idx2], filename)
        print(f"  Created: {filename}")
        print(f"    Contains sequences: {concatenated_seqs[idx1]['source']} and {concatenated_seqs[idx2]['source']}")

if __name__ == "__main__":
    main()