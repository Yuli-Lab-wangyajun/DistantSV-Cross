#!/usr/bin/env python3
"""
Extract BED coordinates from MAF alignments with upstream/downstream extension.
Generate sequence extraction, length comparison, and statistics.
"""

import argparse
import sys
import re

def parse_maf_blocks(maf_file):
    """Parse MAF file into alignment blocks"""
    blocks = []
    current_block = None
    
    with open(maf_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('a'):
                if current_block is not None:
                    blocks.append(current_block)
                current_block = {
                    'header': line,
                    'sequences': []
                }
            elif line.startswith('s') and current_block is not None:
                parts = re.split(r'\s+', line)
                if len(parts) < 7:
                    continue
                    
                seq_info = {
                    'original_line': line,
                    'full_name': parts[1],
                    'start': int(parts[2]),
                    'length': int(parts[3]),
                    'strand': parts[4],
                    'src_length': int(parts[5]),
                    'sequence': parts[6]
                }
                current_block['sequences'].append(seq_info)
        
        if current_block is not None:
            blocks.append(current_block)
    
    return blocks

def extract_species_and_chrom(full_name):
    """Extract species and chromosome name from full sequence name"""
    if '.' in full_name:
        parts = full_name.split('.', 1)
        return parts[0], parts[1]
    return None, full_name

def find_exact_block_for_bed(bed_entry, maf_blocks, ref_species):
    """Find the exact MAF block containing the BED coordinates"""
    chrom = bed_entry['chrom']
    start = bed_entry['start']
    end = bed_entry['end']
    
    for i, block in enumerate(maf_blocks):
        ref_seq_info = None
        
        for seq_info in block['sequences']:
            species, block_chrom = extract_species_and_chrom(seq_info['full_name'])
            if species == ref_species and block_chrom == chrom:
                ref_seq_info = seq_info
                break
        
        if ref_seq_info is None:
            continue
            
        block_start = ref_seq_info['start']
        block_end = block_start + ref_seq_info['length']
        
        if start >= block_start and end <= block_end:
            print(f"Exact match: {chrom}:{start}-{end} -> block index {i} (reference species: {ref_species})")
            return i, ref_seq_info
    
    return None, None

def find_alignment_columns_for_bed(bed_entry, ref_seq_info):
    """
    Find the column range in alignment for BED coordinates
    Returns: (start_col, end_col)
    """
    ref_sequence = ref_seq_info['sequence']
    block_start = ref_seq_info['start']
    start_pos = bed_entry['start']
    end_pos = bed_entry['end']
    bed_length = end_pos - start_pos
    
    # Find start column: track genome position until equal to start_pos
    genome_pos = block_start
    start_col = 0
    
    # Find start column
    while start_col < len(ref_sequence) and genome_pos < start_pos:
        if ref_sequence[start_col] != '-':
            genome_pos += 1
        start_col += 1
    
    # If start position not found
    if start_col >= len(ref_sequence):
        print(f"  Error: Cannot find start position {start_pos} in alignment sequence")
        return None, None
    
    # Calculate end column (based on BED length)
    end_col = start_col + bed_length
    
    # Ensure not exceeding sequence boundary
    if end_col > len(ref_sequence):
        print(f"  Warning: Extraction range exceeds sequence boundary")
        end_col = len(ref_sequence)
    
    return start_col, end_col

def extract_extended_sequences(bed_entry, ref_seq_info, target_seq_info, extend_bp):
    """
    Extract extended sequences
    extend_bp: upstream/downstream extension in base pairs (genome coordinates, excluding gaps)
    """
    ref_sequence = ref_seq_info['sequence']
    target_sequence = target_seq_info['sequence']
    
    # Ensure two sequences have same length
    if len(ref_sequence) != len(target_sequence):
        print(f"  Error: Reference sequence length({len(ref_sequence)}) differs from target sequence length({len(target_sequence)})")
        return None, None, 0, 0, 0, 0
    
    # Find original BED coordinate column range
    start_col, end_col = find_alignment_columns_for_bed(bed_entry, ref_seq_info)
    if start_col is None:
        return None, None, 0, 0, 0, 0
    
    # Upstream extension
    extended_start_col = start_col
    if extend_bp > 0:
        genome_pos = ref_seq_info['start']
        upstream_bp_count = 0
        current_col = start_col - 1
        
        # Traverse left, count non-gap characters
        while current_col >= 0 and upstream_bp_count < extend_bp:
            if ref_sequence[current_col] != '-':
                upstream_bp_count += 1
            current_col -= 1
        
        # Calculate extended start column
        extended_start_col = current_col + 1
        print(f"  Upstream extension: from column {start_col} to {extended_start_col}, extended {upstream_bp_count}bp")
    
    # Downstream extension
    extended_end_col = end_col
    if extend_bp > 0:
        downstream_bp_count = 0
        current_col = end_col
        
        # Traverse right, count non-gap characters
        while current_col < len(ref_sequence) and downstream_bp_count < extend_bp:
            if ref_sequence[current_col] != '-':
                downstream_bp_count += 1
            current_col += 1
        
        # Calculate extended end column
        extended_end_col = current_col
        print(f"  Downstream extension: from column {end_col} to {extended_end_col}, extended {downstream_bp_count}bp")
    
    # Extract extended sequences
    ref_seq = ref_sequence[extended_start_col:extended_end_col]
    target_seq = target_sequence[extended_start_col:extended_end_col]
    
    # Calculate statistics
    ref_total_length = len(ref_seq)  # Alignment length (including gaps)
    target_total_length = len(target_seq)  # Alignment length (including gaps)
    ref_base_count = sum(1 for c in ref_seq if c != '-')  # Genome length (excluding gaps)
    target_base_count = sum(1 for c in target_seq if c != '-')  # Genome length (excluding gaps)
    
    return ref_seq, target_seq, ref_total_length, ref_base_count, target_total_length, target_base_count

def extract_sequences_precise(bed_entry, ref_seq_info, target_seq_info):
    """
    Extract BED coordinate sequences precisely (no extension)
    """
    ref_sequence = ref_seq_info['sequence']
    target_sequence = target_seq_info['sequence']
    
    if len(ref_sequence) != len(target_sequence):
        print(f"  Error: Reference sequence length({len(ref_sequence)}) differs from target sequence length({len(target_sequence)})")
        return None, None, 0, 0, 0, 0
    
    # Find alignment column range
    start_col, end_col = find_alignment_columns_for_bed(bed_entry, ref_seq_info)
    if start_col is None:
        return None, None, 0, 0, 0, 0
    
    # Extract sequences
    ref_seq = ref_sequence[start_col:end_col]
    target_seq = target_sequence[start_col:end_col]
    
    # Calculate statistics
    ref_total_length = len(ref_seq)
    target_total_length = len(target_seq)
    ref_base_count = sum(1 for c in ref_seq if c != '-')
    target_base_count = sum(1 for c in target_seq if c != '-')
    
    return ref_seq, target_seq, ref_total_length, ref_base_count, target_total_length, target_base_count

def main():
    parser = argparse.ArgumentParser(description='Extract BED coordinates from MAF alignments with extension')
    parser.add_argument('--bed', required=True, help='Input BED file path')
    parser.add_argument('--maf', required=True, help='Input MAF file path')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--ref-species', type=str, required=True,
                       help='Reference species, e.g., "papo"')
    parser.add_argument('--target-species', type=str, required=True,
                       help='Target species, e.g., "tbai"')
    parser.add_argument('--extend', type=int, default=0,
                       help='Upstream/downstream extension in base pairs (genome coordinates, excluding gaps), default 0')
    parser.add_argument('--threshold', type=float, default=0.5,
                       help='Threshold for large difference between actual check length and detected length (ratio of detected length), default 0.5')
    
    args = parser.parse_args()
    
    print("Parsing BED file...")
    bed_entries = parse_bed_file(args.bed)
    print(f"Found {len(bed_entries)} BED entries")
    
    print("Parsing MAF file...")
    maf_blocks = parse_maf_blocks(args.maf)
    print(f"Found {len(maf_blocks)} MAF alignment blocks")
    
    print(f"Reference species: {args.ref_species}, Target species: {args.target_species}")
    print(f"Extension parameter: {args.extend}bp upstream/downstream" if args.extend > 0 else "Extension parameter: no extension")
    print(f"Threshold parameter: {args.threshold}")
    print("Matching BED coordinates and extracting sequences...")
    
    output_lines = []
    output_lines.append("# BED coordinate extraction results")
    if args.extend > 0:
        output_lines.append(f"# Extended {args.extend}bp upstream/downstream")
    output_lines.append("#" + "="*80)
    
    # Statistics
    total_entries = len(bed_entries)
    found_entries = 0
    not_found_entries = []
    length_comparison_lines = ["Detected_length\tActual_check_length\tLarge_difference"]  # Length comparison file header
    large_diff_count = 0  # Count of entries with large difference
    
    target_lengths_lines = []
    extend_double = 2 * args.extend
    
    for i, bed_entry in enumerate(bed_entries):
        print(f"\nProcessing BED entry {i+1}/{len(bed_entries)}: {bed_entry['chrom']}:{bed_entry['start']}-{bed_entry['end']}")
        bed_length = bed_entry['end'] - bed_entry['start']
        print(f"  BED length: {bed_length}bp")
        
        # Find MAF block containing BED coordinates
        block_idx, ref_seq_info = find_exact_block_for_bed(bed_entry, maf_blocks, args.ref_species)
        
        if block_idx is None or ref_seq_info is None:
            print(f"  Not found: No MAF block contains BED coordinates")
            not_found_entries.append(bed_entry['original_line'])
            # For not found entries, record detected length and NA
            length_comparison_lines.append(f"{bed_length}\tNA\tNA")
            continue
        
        # Find target species sequence
        target_seq_info = None
        block = maf_blocks[block_idx]
        for seq_info in block['sequences']:
            species, _ = extract_species_and_chrom(seq_info['full_name'])
            if species == args.target_species:
                target_seq_info = seq_info
                break
        
        if target_seq_info is None:
            print(f"  Warning: Target species {args.target_species} not found in alignment block {block_idx}")
            not_found_entries.append(bed_entry['original_line'])
            # For entries without target species, record detected length and NA
            length_comparison_lines.append(f"{bed_length}\tNA\tNA")
            continue
        
        if args.extend > 0:
            # Extract extended sequences
            ref_seq, target_seq, ref_total_len, ref_base_count, target_total_len, target_base_count = extract_extended_sequences(
                bed_entry, ref_seq_info, target_seq_info, args.extend
            )
            
            if ref_seq is None:
                print(f"  Warning: Cannot extract extended sequences")
                not_found_entries.append(bed_entry['original_line'])
                length_comparison_lines[-1] = f"{bed_length}\tNA\tNA"  # Update to NA
                large_diff_count -= 1  # Subtract previously counted large difference
                continue
            
            # Calculate actual check length
            # Actual check length = Reference alignment length - Target genome length
            actual_check_length = ref_total_len - target_base_count
            print(f"  Actual check length: {ref_total_len} - {target_base_count} = {actual_check_length}")
            
            # Check if difference is large
            diff_ratio = abs(actual_check_length - bed_length) / bed_length if bed_length > 0 else 0
            diff_large = diff_ratio > args.threshold
            diff_mark = "Yes" if diff_large else "No"
            
            if diff_large:
                large_diff_count += 1
                print(f"  Warning: Large difference between detected length({bed_length}) and actual check length({actual_check_length}), ratio={diff_ratio:.2f}")
            
            # Record length comparison information
            length_comparison_lines.append(f"{bed_length}\t{actual_check_length}\t{diff_mark}")
            
            target_lengths_lines.append(f"{extend_double}\t{target_base_count}")
            # ---------------------------------
            
            # Add to output
            output_lines.append("")
            output_lines.append(f"# BED entry: {bed_entry['original_line']}")
            output_lines.append(f"# Alignment block index: {block_idx}")
            output_lines.append(f"# Extension: original {bed_length} columns, extended to {len(ref_seq)} columns")
            output_lines.append(f"# Reference species sequence ({args.ref_species}, genome coordinates {bed_entry['start']}-{bed_entry['end']}, extended {args.extend}bp):")
            output_lines.append(ref_seq)
            output_lines.append(f"# Target species sequence ({args.target_species}, same alignment columns):")
            output_lines.append(target_seq)
            output_lines.append(f"# Sequence statistics: ")
            output_lines.append(f"#   Reference species: alignment length={ref_total_len} (with gaps), genome length={ref_base_count} (without gaps)")
            output_lines.append(f"#   Target species: alignment length={target_total_len} (with gaps), genome length={target_base_count} (without gaps)")
            
            print(f"  Extended extraction successful: reference {ref_total_len} chars (with gaps), {ref_base_count}bp (without gaps)")
            print(f"                                target {target_total_len} chars (with gaps), {target_base_count}bp (without gaps)")
        else:
            # Extract original sequences (no extension)
            ref_seq, target_seq, ref_total_len, ref_base_count, target_total_len, target_base_count = extract_sequences_precise(
                bed_entry, ref_seq_info, target_seq_info
            )
            
            if ref_seq is None:
                print(f"  Warning: Cannot extract sequences")
                not_found_entries.append(bed_entry['original_line'])
                length_comparison_lines[-1] = f"{bed_length}\tNA\tNA"  # Update to NA
                large_diff_count -= 1  # Subtract previously counted large difference
                continue
            
            # Calculate actual check length
            # Actual check length = Reference alignment length - Target genome length
            actual_check_length = ref_total_len - target_base_count
            print(f"  Actual check length: {ref_total_len} - {target_base_count} = {actual_check_length}")
            
            # Check if difference is large
            diff_ratio = abs(actual_check_length - bed_length) / bed_length if bed_length > 0 else 0
            diff_large = diff_ratio > args.threshold
            diff_mark = "Yes" if diff_large else "No"
            
            if diff_large:
                large_diff_count += 1
                print(f"  Warning: Large difference between detected length({bed_length}) and actual check length({actual_check_length}), ratio={diff_ratio:.2f}")
            
            # Record length comparison information
            length_comparison_lines.append(f"{bed_length}\t{actual_check_length}\t{diff_mark}")
            
            target_lengths_lines.append(f"{extend_double}\t{target_base_count}")
            # -----------------------------------------------------------
            
            # Add to output
            output_lines.append("")
            output_lines.append(f"# BED entry: {bed_entry['original_line']}")
            output_lines.append(f"# Alignment block index: {block_idx}")
            output_lines.append(f"# Reference species sequence ({args.ref_species}, genome coordinates {bed_entry['start']}-{bed_entry['end']}):")
            output_lines.append(ref_seq)
            output_lines.append(f"# Target species sequence ({args.target_species}, same alignment columns):")
            output_lines.append(target_seq)
            output_lines.append(f"# Sequence statistics: ")
            output_lines.append(f"#   Reference species: alignment length={ref_total_len} (with gaps), genome length={ref_base_count} (without gaps)")
            output_lines.append(f"#   Target species: alignment length={target_total_len} (with gaps), genome length={target_base_count} (without gaps)")
            
            print(f"  Extraction successful: reference {ref_total_len} chars, {ref_base_count}bp (without gaps)")
            print(f"                        target {target_total_len} chars, {target_base_count}bp (without gaps)")
        
        found_entries += 1
    
    # Write output file
    with open(args.output, 'w') as f:
        f.write('\n'.join(output_lines))
    
    # Write length comparison file
    length_comparison_file = args.output + '.length_comparison.txt'
    with open(length_comparison_file, 'w') as f:
        f.write('\n'.join(length_comparison_lines))
    
    # Write statistics file
    stats_file = args.output + '.statistics.txt'
    not_found_count = len(not_found_entries)
    
    # Calculate indeterminate entries (actual check length is NA)
    indeterminate_count = 0
    for line in length_comparison_lines[1:]:
        parts = line.split('\t')
        if len(parts) >= 2 and parts[1] == 'NA':
            indeterminate_count += 1
    
    stats_lines = [
        "# Statistics",
        "#" + "="*40,
        f"Total entries: {total_entries}",
        f"Entries with matching alignment blocks: {found_entries}",
        f"Entries without matching alignment blocks: {not_found_count}",
        f"Indeterminate entries: {indeterminate_count}",
        f"Entries with large difference: {large_diff_count} (threshold: {args.threshold})",
    ]
    
    if total_entries > 0:
        stats_lines.append(f"Large difference entries / Total entries: {large_diff_count/total_entries:.2%}")
    else:
        stats_lines.append(f"Large difference entries / Total entries: 0.00%")
    
    if found_entries > 0:
        stats_lines.append(f"Large difference entries / Matched entries: {large_diff_count/found_entries:.2%}")
    else:
        stats_lines.append(f"Large difference entries / Matched entries: 0.00%")
    
    with open(stats_file, 'w') as f:
        f.write('\n'.join(stats_lines))
    
    # Handle not found entries
    if not_found_entries:
        not_found_file = args.output + '.not_found'
        print(f"\nWarning: {not_found_count} BED entries have no matching MAF alignment blocks")
        print(f"Not found entries written to: {not_found_file}")
        
        with open(not_found_file, 'w') as f:
            f.write('\n'.join(not_found_entries))
    
    target_lengths_file = args.output + '.target_lengths.txt'
    with open(target_lengths_file, 'w') as f:
        f.write("# extend_double\ttarget_genome_length\n")
        f.write('\n'.join(target_lengths_lines))
    
    print(f"\nSuccessfully extracted sequences for {found_entries} BED entries")
    print(f"Length comparison file generated: {length_comparison_file}")
    print(f"Statistics file generated: {stats_file}")
    print(f"Target genome lengths file generated: {target_lengths_file}")
    print("Processing completed!")

def parse_bed_file(bed_file):
    """Parse BED file"""
    bed_entries = []
    with open(bed_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            parts = line.split('\t')
            if len(parts) < 3:
                continue
                
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            bed_entries.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'original_line': line
            })
    
    return bed_entries

if __name__ == '__main__':
    main()