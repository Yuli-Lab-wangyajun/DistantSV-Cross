#!/usr/bin/env python3
"""
Loosely match the position of BED coordinates in MAF files, specify the species to check,
match if BED start or end is within the alignment block, output all matching blocks
"""

import argparse
import sys
import re

def parse_maf_blocks(maf_file):
    """Parse MAF file into a list of alignment blocks, support tab and space separation"""
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
                # Parse sequence line: s species.chrom start length strand src_length sequence
                # Support tab and space separation
                parts = re.split(r'\s+', line)
                if len(parts) < 7:
                    continue
                    
                seq_info = {
                    'original_line': line,
                    'full_name': parts[1],      # Full name, e.g.: tbai.QLTV01002257.1
                    'start': int(parts[2]),     # Sequence start position
                    'length': int(parts[3]),    # Sequence length
                    'strand': parts[4],         # Strand direction
                    'src_length': int(parts[5]), # Total length of source sequence
                    'sequence': parts[6]        # Alignment sequence
                }
                current_block['sequences'].append(seq_info)
        
        if current_block is not None:
            blocks.append(current_block)
    
    return blocks

def extract_species_and_chrom(full_name):
    """
    Extract species and chromosome name from full sequence name
    Example: tbai.QLTV01002257.1 -> (tbai, QLTV01002257.1)
    """
    if '.' in full_name:
        parts = full_name.split('.', 1)  # Split only at the first dot
        return parts[0], parts[1]
    return None, full_name  # If no dot, return None and full name

def find_overlapping_blocks_for_bed(bed_entry, maf_blocks, target_species=None):
    """
    Find MAF alignment blocks overlapping with BED coordinates
    
    Parameters:
    - bed_entry: BED entry information
    - maf_blocks: List of MAF alignment blocks
    - target_species: Specify species to check, e.g. 'papo' or 'tbai'
    
    Matching conditions:
    - BED start position is within MAF block range, or
    - BED end position is within MAF block range, or
    - BED completely contains MAF block, or
    - MAF block completely contains BED
    """
    chrom = bed_entry['chrom']
    start = bed_entry['start']
    end = bed_entry['end']
    
    matching_blocks = []
    
    for i, block in enumerate(maf_blocks):
        for seq_info in block['sequences']:
            # Extract species and chromosome name
            species, block_chrom = extract_species_and_chrom(seq_info['full_name'])
            
            # If species is specified, only check sequences of that species
            if target_species and species != target_species:
                continue
                
            # Check if chromosome name matches
            if block_chrom == chrom:
                block_start = seq_info['start']
                block_end = block_start + seq_info['length']
                
                # Loose matching condition: match if BED start or end is within MAF block, or they overlap
                if ((start >= block_start and start <= block_end) or  # BED start within MAF block
                    (end >= block_start and end <= block_end) or      # BED end within MAF block
                    (start <= block_start and end >= block_end) or    # BED completely contains MAF block
                    (start >= block_start and end <= block_end)):     # MAF block completely contains BED
                    
                    matching_blocks.append(i)
                    print(f"Match found: {chrom}:{start}-{end} -> Block index {i} (Species: {species})")
                    print(f"  BED range: {start}-{end}")
                    print(f"  MAF block range: {block_start}-{block_end}")
                    break  # One matching sequence per block is sufficient
    
    return matching_blocks

def find_blocks_in_range(bed_entry, maf_blocks, search_range=1000, target_species=None):
    """
    If not found in loose matching, further relax conditions to search nearby
    """
    chrom = bed_entry['chrom']
    start = bed_entry['start']
    end = bed_entry['end']
    
    matching_blocks = []
    
    for i, block in enumerate(maf_blocks):
        for seq_info in block['sequences']:
            species, block_chrom = extract_species_and_chrom(seq_info['full_name'])
            
            # If species is specified, only check sequences of that species
            if target_species and species != target_species:
                continue
                
            if block_chrom == chrom:
                block_start = seq_info['start']
                block_end = block_start + seq_info['length']
                
                # Further relax matching condition: match if overlapping within specified range
                expanded_start = block_start - search_range
                expanded_end = block_end + search_range
                
                if ((start >= expanded_start and start <= expanded_end) or
                    (end >= expanded_start and end <= expanded_end) or
                    (start <= expanded_start and end >= expanded_end) or
                    (start >= expanded_start and end <= expanded_end)):
                    
                    matching_blocks.append(i)
                    print(f"Range match: {chrom}:{start}-{end} -> Block index {i} (Species: {species})")
                    print(f"  BED range: {start}-{end}")
                    print(f"  MAF block range: {block_start}-{block_end}")
                    print(f"  Expanded range: {expanded_start}-{expanded_end}")
                    break
    
    return matching_blocks

def analyze_maf_coverage(maf_blocks, target_chrom, target_species=None):
    """
    Analyze coverage of specific chromosome in MAF file, for debugging
    """
    coverage = []
    for i, block in enumerate(maf_blocks):
        for seq_info in block['sequences']:
            species, chrom = extract_species_and_chrom(seq_info['full_name'])
            
            # If species is specified, only check sequences of that species
            if target_species and species != target_species:
                continue
                
            if chrom == target_chrom:
                coverage.append({
                    'block_index': i,
                    'species': species,
                    'chrom': chrom,
                    'start': seq_info['start'],
                    'end': seq_info['start'] + seq_info['length'],
                    'length': seq_info['length']
                })
    
    return coverage

def main():
    parser = argparse.ArgumentParser(description='Loosely match BED coordinates position in MAF files')
    parser.add_argument('--bed', required=True, help='Input BED file path')
    parser.add_argument('--maf', required=True, help='Input MAF file path')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--context', type=int, default=1, help='Number of adjacent alignment blocks to extract')
    parser.add_argument('--search-range', type=int, default=0, 
                       help='Search range, 0 for loose matching, >0 for search within specified range')
    parser.add_argument('--species', type=str, default=None,
                       help='Specify species to check, e.g. "papo" or "tbai". Use all species if not specified')
    
    args = parser.parse_args()
    
    print("Parsing BED file...")
    bed_entries = parse_bed_file(args.bed)
    print(f"Found {len(bed_entries)} BED entries")
    
    print("Parsing MAF file...")
    maf_blocks = parse_maf_blocks(args.maf)
    print(f"Found {len(maf_blocks)} MAF alignment blocks")
    
    # Analyze chromosome coverage of first BED entry (for debugging)
    if bed_entries:
        first_chrom = bed_entries[0]['chrom']
        coverage = analyze_maf_coverage(maf_blocks, first_chrom, args.species)
        species_info = f" (Species: {args.species})" if args.species else " (All species)"
        print(f"Coverage of chromosome {first_chrom}{species_info} in MAF: {len(coverage)} blocks")
        if coverage:
            # Sort by start position
            coverage.sort(key=lambda x: x['start'])
            print(f"  Coordinate range: {coverage[0]['start']}-{coverage[-1]['end']}")
            print(f"  Total coverage length: {coverage[-1]['end'] - coverage[0]['start']}")
            if coverage:
                print(f"  Example species: {coverage[0]['species']}")
                
                # Show details of first 5 blocks
                print(f"  First 5 blocks:")
                for i, cov in enumerate(coverage[:5]):
                    print(f"    Block {cov['block_index']}: {cov['start']}-{cov['end']} (Length: {cov['length']})")
    
    print("Matching BED coordinates with MAF alignment blocks...")
    not_found_entries = []
    output_lines = []
    
    for i, bed_entry in enumerate(bed_entries):
        if args.search_range == 0:
            matching_blocks = find_overlapping_blocks_for_bed(bed_entry, maf_blocks, args.species)
        else:
            matching_blocks = find_blocks_in_range(bed_entry, maf_blocks, args.search_range, args.species)
        
        if not matching_blocks:
            species_info = f" (Species: {args.species})" if args.species else ""
            print(f"No match found: {bed_entry['chrom']}:{bed_entry['start']}-{bed_entry['end']}{species_info}")
            not_found_entries.append(bed_entry['original_line'])
            continue
        
        # Deduplicate matching block indices
        matching_blocks = list(set(matching_blocks))
        matching_blocks.sort()
        
        # Extract matching blocks and their context
        context_blocks = extract_blocks_with_context(matching_blocks, maf_blocks, args.context)
        
        # Build output - add comment line indicating number of matched blocks
        output_lines.append(bed_entry['original_line'])
        output_lines.append(f"# This BED entry matched {len(matching_blocks)} alignment blocks")
        
        for j, block in enumerate(context_blocks):
            output_lines.append(f"# Alignment block {j+1}/{len(context_blocks)}")
            output_lines.append(block['header'])
            for seq_info in block['sequences']:
                output_lines.append(seq_info['original_line'])
        
        output_lines.append('')
        output_lines.append('')
    
    # Write to output file
    with open(args.output, 'w') as f:
        f.write('\n'.join(output_lines))
    
    # Process entries not found
    if not_found_entries:
        not_found_file = args.output + '.not_found'
        species_info = f" (Species: {args.species})" if args.species else ""
        print(f"Warning: {len(not_found_entries)} BED entries did not find corresponding MAF alignment blocks{species_info}")
        print(f"Entries not found have been written to: {not_found_file}")
        
        with open(not_found_file, 'w') as f:
            f.write('\n'.join(not_found_entries))
        
        # Suggest adjusting search range
        print("Suggestion: If matches should exist but were not found, try expanding search range with --search-range parameter")
        if args.species:
            print(f"      Or check if --species parameter is correct (Current: {args.species})")
        
        # Provide more debugging information
        if bed_entries:
            first_not_found = not_found_entries[0].split('\t')
            if len(first_not_found) >= 3:
                chrom = first_not_found[0]
                start = int(first_not_found[1])
                end = int(first_not_found[2])
                print(f"Details of example entry not found:")
                print(f"  BED coordinate: {chrom}:{start}-{end}")
                
                # Find all blocks for this chromosome
                coverage = analyze_maf_coverage(maf_blocks, chrom, args.species)
                if coverage:
                    print(f"  Coverage range of this chromosome in MAF: {coverage[0]['start']}-{coverage[-1]['end']}")
                    # Find closest block
                    closest_block = min(coverage, key=lambda x: min(abs(x['start'] - start), abs(x['end'] - end)))
                    print(f"  Closest block: Index {closest_block['block_index']}, Range {closest_block['start']}-{closest_block['end']}")
                    print(f"  Distance: {min(abs(closest_block['start'] - start), abs(closest_block['end'] - end))} bp")
    else:
        species_info = f" (Species: {args.species})" if args.species else ""
        print(f"All BED entries successfully found corresponding MAF alignment blocks{species_info}")
    
    print("Processing completed!")

def extract_blocks_with_context(block_indices, maf_blocks, context=1):
    """Extract specified blocks and their adjacent context blocks, deduplicate"""
    if not block_indices:
        return []
    
    result_blocks = []
    all_indices = set()
    
    for idx in block_indices:
        for i in range(max(0, idx - context), min(len(maf_blocks), idx + context + 1)):
            all_indices.add(i)
    
    sorted_indices = sorted(all_indices)
    return [maf_blocks[i] for i in sorted_indices]

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