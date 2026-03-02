#!/usr/bin/env python3
"""
Calculate conservation metrics for sequences flanking missing variants
Input: MAF multiple sequence alignment file and BED format variant file
Output: BED file with added conservation metrics
"""

import argparse
import sys
from typing import List, Dict, Tuple, Optional
import pandas as pd
from collections import Counter

class MafBlock:
    """MAF alignment block class"""
    def __init__(self, score: float = 0.0):
        self.score = score
        self.sequences = []  # Each sequence tuple: (species, chrom, start, length, strand, seq)
    
    def add_sequence(self, line_parts: List[str]):
        """Add sequence to alignment block"""
        # Parse sequence line: s species.chrom start length strand total_length sequence
        name_parts = line_parts[1].split('.', 1)
        species = name_parts[0]
        chrom = name_parts[1] if len(name_parts) > 1 else name_parts[0]
        start = int(line_parts[2])
        length = int(line_parts[3])
        strand = line_parts[4]
        seq = line_parts[6]
        self.sequences.append((species, chrom, start, length, strand, seq))
    
    def get_alignment_matrix(self) -> List[List[str]]:
        """Get alignment matrix (one list per column)"""
        if not self.sequences:
            return []
        
        seq_length = len(self.sequences[0][5])
        matrix = []
        for i in range(seq_length):
            column = [seq[5][i].upper() for seq in self.sequences]
            matrix.append(column)
        return matrix
    
    def find_reference_position(self, ref_species: str, ref_pos: int) -> Optional[int]:
        """Find column index in alignment matrix for specific reference sequence position"""
        for i, (species, chrom, start, length, strand, seq) in enumerate(self.sequences):
            if species == ref_species:
                if start <= ref_pos < start + length:
                    # Calculate column position in alignment
                    pos_in_seq = ref_pos - start
                    # Need to account for gaps in sequence, find corresponding alignment column
                    col_idx = 0
                    seq_idx = 0
                    while seq_idx <= pos_in_seq and col_idx < len(seq):
                        if seq[col_idx] != '-':  # Non-gap position
                            seq_idx += 1
                        if seq_idx > pos_in_seq:
                            break
                        col_idx += 1
                    return col_idx if col_idx < len(seq) else None
        return None

def read_maf_file(maf_file: str) -> List[MafBlock]:
    """Read MAF file"""
    blocks = []
    current_block = None
    
    with open(maf_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('a score='):
                if current_block is not None:
                    blocks.append(current_block)
                score = float(line.split('=')[1].split()[0])
                current_block = MafBlock(score)
            elif line.startswith('s ') and current_block is not None:
                parts = line.split()
                if len(parts) >= 7:
                    current_block.add_sequence(parts)
    
    if current_block is not None:
        blocks.append(current_block)
    
    return blocks

def calculate_column_identity(column: List[str]) -> float:
    """
    Calculate identity percentage for a single column
    Denominator is always total number of sequences, numerator is count of identical bases (excluding gaps)
    """
    if not column:
        return 0.0
    
    total_sequences = len(column)
    
    # Filter out gap characters, only consider non-gap bases
    non_gap_bases = [base for base in column if base != '-']
    
    # If no non-gap bases, identity is 0
    if not non_gap_bases:
        return 0.0
    
    # Find most common non-gap base
    base_counts = Counter(non_gap_bases)
    most_common_count = base_counts.most_common(1)[0][1]
    
    # Identity = count of most common base / total sequences
    return most_common_count / total_sequences

def calculate_average_identity(columns: List[List[str]]) -> float:
    """Calculate average identity across multiple columns"""
    if not columns:
        return 0.0
    
    total_identity = sum(calculate_column_identity(col) for col in columns)
    return total_identity / len(columns)

def get_conservation_for_intra_block_deletion(block: MafBlock, ref_species: str, 
                                            del_start: int, del_end: int) -> Tuple[float, float, float, int, int]:
    """
    Calculate conservation for deletion within a block
    Use entire alignment block (except deletion region) for calculation
    Returns: (upstream conservation, downstream conservation, total conservation, upstream length, downstream length)
    """
    # Find column indices in alignment for deletion positions in reference sequence
    del_start_col = block.find_reference_position(ref_species, del_start)
    del_end_col = block.find_reference_position(ref_species, del_end - 1)  # BED format is half-open interval
    
    if del_start_col is None or del_end_col is None:
        return 0.0, 0.0, 0.0, 0, 0
    
    # Get alignment matrix
    matrix = block.get_alignment_matrix()
    total_columns = len(matrix)
    
    # Calculate upstream conservation (all columns before deletion)
    upstream_start = 0
    upstream_end = del_start_col
    upstream_columns = matrix[upstream_start:upstream_end]
    upstream_conservation = calculate_average_identity(upstream_columns) if upstream_columns else 0.0
    upstream_length = len(upstream_columns)
    
    # Calculate downstream conservation (all columns after deletion)
    downstream_start = del_end_col + 1  # Skip deletion region
    downstream_end = total_columns
    downstream_columns = matrix[downstream_start:downstream_end]
    downstream_conservation = calculate_average_identity(downstream_columns) if downstream_columns else 0.0
    downstream_length = len(downstream_columns)
    
    # Calculate total conservation (average of upstream and downstream)
    total_conservation = (upstream_conservation + downstream_conservation) / 2
    
    return upstream_conservation, downstream_conservation, total_conservation, upstream_length, downstream_length

def get_conservation_for_inter_block_deletion(blocks: List[MafBlock]) -> Tuple[float, float, float, int, int]:
    """
    Calculate conservation for deletion between blocks
    For multiple blocks, only use first and last block
    Returns: (upstream block conservation, downstream block conservation, total conservation, upstream length, downstream length)
    """
    if not blocks:
        return 0.0, 0.0, 0.0, 0, 0
    
    # Only use first and last block
    first_block = blocks[0]
    last_block = blocks[-1]
    
    # Calculate conservation for upstream block (first block)
    matrix1 = first_block.get_alignment_matrix()
    upstream_conservation = calculate_average_identity(matrix1) if matrix1 else 0.0
    upstream_length = len(matrix1) if matrix1 else 0
    
    # Calculate conservation for downstream block (last block)
    matrix2 = last_block.get_alignment_matrix()
    downstream_conservation = calculate_average_identity(matrix2) if matrix2 else 0.0
    downstream_length = len(matrix2) if matrix2 else 0
    
    # Calculate total conservation (average of two blocks)
    total_conservation = (upstream_conservation + downstream_conservation) / 2
    
    return upstream_conservation, downstream_conservation, total_conservation, upstream_length, downstream_length

def main():
    parser = argparse.ArgumentParser(description='Calculate conservation metrics for sequences flanking missing variants')
    parser.add_argument('-m', '--maf', required=True, help='Input MAF file')
    parser.add_argument('-b', '--bed', required=True, help='Input BED file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')
    parser.add_argument('-r', '--ref-species', default='mu_mu', help='Reference species name (default: mu_mu)')
    
    args = parser.parse_args()
    
    # Read MAF file
    print("Reading MAF file...", file=sys.stderr)
    maf_blocks = read_maf_file(args.maf)
    print(f"Read {len(maf_blocks)} alignment blocks", file=sys.stderr)
    
    # Read BED file - directly use known column names
    print("Reading BED file...", file=sys.stderr)
    
    # Define BED file column names
    bed_columns = [
        'Ref_Chrom', 'Ref_Start', 'Ref_End', 'Ref_Strand',
        'Target_Chrom', 'Target_Start', 'Target_End', 'Target_Strand',
        'Original_tra_Chrom', 'Original_tra_Start', 'Original_tra_End', 'Original_tra_Strand',
        'SV_Type', 'SV_Size', 'SV_Source', 'Block_Indexes', 'Additional_Info'
    ]
    
    # Read BED file, skip comment lines, use custom column names, and preserve "NA" strings
    bed_df = pd.read_csv(
        args.bed, 
        sep='\t', 
        comment='#', 
        header=None, 
        names=bed_columns,
        keep_default_na=False,  # Don't convert "NA" strings to NaN
        na_values=[]  # Don't treat any values as NaN
    )
    
    # Print BED file column names for debugging
    print(f"BED file column names: {list(bed_df.columns)}", file=sys.stderr)
    print(f"BED file row count: {len(bed_df)}", file=sys.stderr)
    
    # Add new columns for storing conservation results and length information
    bed_df['SV_Upstream_Conservation'] = 0.0
    bed_df['SV_Downstream_Conservation'] = 0.0
    bed_df['SV_Total_Conservation'] = 0.0
    bed_df['SV_Upstream_Length'] = 0
    bed_df['SV_Downstream_Length'] = 0
    
    # Process each variant
    for idx, row in bed_df.iterrows():
        block_indexes = [int(idx) for idx in str(row['Block_Indexes']).split(',')]
        sv_source = row['SV_Source']
        
        try:
            if sv_source == 'intra-block':
                # Deletion within a single block
                if len(block_indexes) >= 1:
                    block_idx = block_indexes[0]
                    if block_idx < len(maf_blocks):
                        block = maf_blocks[block_idx]
                        upstream_cons, downstream_cons, total_cons, upstream_len, downstream_len = get_conservation_for_intra_block_deletion(
                            block, args.ref_species, int(row['Ref_Start']), int(row['Ref_End'])
                        )
                        bed_df.at[idx, 'SV_Upstream_Conservation'] = upstream_cons
                        bed_df.at[idx, 'SV_Downstream_Conservation'] = downstream_cons
                        bed_df.at[idx, 'SV_Total_Conservation'] = total_cons
                        bed_df.at[idx, 'SV_Upstream_Length'] = upstream_len
                        bed_df.at[idx, 'SV_Downstream_Length'] = downstream_len
                    else:
                        print(f"Warning: Block index {block_idx} out of range", file=sys.stderr)
                else:
                    print(f"Warning: intra-block requires at least 1 block index, but only {len(block_indexes)} found", file=sys.stderr)
            
            elif sv_source == 'inter-block':
                # Deletion between two or more blocks
                if len(block_indexes) >= 2:
                    # Get all relevant blocks
                    relevant_blocks = [maf_blocks[i] for i in block_indexes if i < len(maf_blocks)]
                    
                    if len(relevant_blocks) >= 2:
                        upstream_cons, downstream_cons, total_cons, upstream_len, downstream_len = get_conservation_for_inter_block_deletion(
                            relevant_blocks
                        )
                        bed_df.at[idx, 'SV_Upstream_Conservation'] = upstream_cons
                        bed_df.at[idx, 'SV_Downstream_Conservation'] = downstream_cons
                        bed_df.at[idx, 'SV_Total_Conservation'] = total_cons
                        bed_df.at[idx, 'SV_Upstream_Length'] = upstream_len
                        bed_df.at[idx, 'SV_Downstream_Length'] = downstream_len
                    else:
                        print(f"Warning: Only {len(relevant_blocks)} valid blocks, need at least 2", file=sys.stderr)
                else:
                    print(f"Warning: inter-block requires at least 2 block indexes, but only {len(block_indexes)} found", file=sys.stderr)
            
            else:
                print(f"Warning: Unknown SV_Source: {sv_source}", file=sys.stderr)
                
        except Exception as e:
            print(f"Error processing row {idx}: {e}", file=sys.stderr)
    
    # Save results
    print(f"Saving results to {args.output}", file=sys.stderr)
    bed_df.to_csv(args.output, sep='\t', index=False)
    print("Completed!", file=sys.stderr)

if __name__ == "__main__":
    main()