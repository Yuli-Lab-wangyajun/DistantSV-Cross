#!/usr/bin/env python3
import argparse
import sys
import os
import re
from collections import defaultdict

def parse_maf_blocks(maf_file):
    """Parse MAF file, return list of blocks, preserving sequence information"""
    blocks = []
    current_block = {}
    
    with open(maf_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Handle empty lines (used to separate blocks in Cactus MAF)
            if not line:
                if current_block:
                    blocks.append(current_block)
                    current_block = {}
                continue
                
            # Handle alignment block start line (a line)
            if line.startswith('a'):
                if current_block:
                    blocks.append(current_block)
                
                # Handle Cactus format (may not have score=)
                if 'score=' in line:
                    score = float(line.split('score=')[1].split()[0])
                else:
                    score = 0.0  # Default value
                
                current_block = {'score': score, 'species': {}}
            
            # Handle sequence line (s line)
            elif line.startswith('s'):
                parts = line.split()
                if len(parts) < 7:  # Ensure enough fields
                    continue
                    
                full_name = parts[1]
                
                # More flexible sequence name parsing
                if '.' in full_name:
                    # Take part before first dot as species name
                    species = full_name.split('.')[0]
                    # Remaining part as chromosome name
                    chrom = '.'.join(full_name.split('.')[1:])
                else:
                    species = full_name
                    chrom = full_name
                
                # Extract sequence information
                try:
                    start = int(parts[2])
                    size = int(parts[3])
                    strand = parts[4]
                    chrom_size = int(parts[5])
                    seq = parts[6]
                except (ValueError, IndexError):
                    continue  # Skip malformed lines
                
                # Calculate actual positions (considering strand direction)
                if strand == '-':
                    actual_start = chrom_size - (start + size)
                    actual_end = chrom_size - start
                else:
                    actual_start = start
                    actual_end = start + size
                
                # Store sequence information
                current_block['species'][species] = {
                    'chrom': chrom,
                    'start': start,
                    'end': start + size,
                    'size': size,
                    'strand': strand,
                    'chrom_size': chrom_size,
                    'seq': seq,
                    'actual_start': actual_start,
                    'actual_end': actual_end
                }
        
        # Process last block
        if current_block:
            blocks.append(current_block)
    
    return blocks

def calculate_gap(prev, curr, species):
    if prev['chrom'] != curr['chrom'] or prev['strand'] != curr['strand']:
        return None

    # Currently MAF coordinates increase linearly, gap calculation same for both strands
    return curr['start'] - prev['end']

def calc_genomic_interval(start, size, strand, seq, start_pos, end_pos, chrom_size):
    """
    Calculate genomic coordinate interval for a continuous region
    """
    # Count valid bases before the region
    valid_before = 0
    for i in range(start_pos):
        if seq[i] != '-':
            valid_before += 1
    
    region_length = end_pos - start_pos
    
    if strand == '+':
        start_coord = start + valid_before
        end_coord = start_coord + region_length
        return start_coord, end_coord
    else:
        # Correct position calculation on negative strand
        valid_in_region = 0
        for i in range(start_pos, end_pos):
            if seq[i] != '-':
                valid_in_region += 1
        
        genome_start = start + valid_before
        start_coord = chrom_size - (genome_start + valid_in_region)
        end_coord = chrom_size - genome_start
        
        # Ensure coordinate order
        if start_coord > end_coord:
            start_coord, end_coord = end_coord, start_coord
        return start_coord, end_coord

def detect_block_svs(blocks, target, ref_species, non_target_species, min_size=10, max_gap_ratio=0.2):
    """Intra-block SV detection: precisely detect target species-specific insertions and deletions"""
    sv_events = []
    
    for block_idx, block in enumerate(blocks):
        species_data = block['species']
        
        # Check if target species and reference species exist
        if target not in species_data or ref_species not in species_data:
            continue
            
        # Get specified non-target species (only those present in current block)
        common_species = set(non_target_species) & set(species_data.keys())
        if not common_species:
            continue
            
        # Get sequence length (should be same for all species)
        seq_len = len(species_data[target]['seq'])
        
        # ================= 1. Detect target species-specific insertions =================
        # Target species has base, while all other species have gaps
        pos = 0
        while pos < seq_len:
            # Check if target species has base at this position
            target_has_base = species_data[target]['seq'][pos] != '-'
            
            # Check if all non-target species have gaps at this position
            non_target_gap = all(
                species_data[sp]['seq'][pos] == '-' 
                for sp in common_species
            )
            
            # Found insertion start point
            if target_has_base and non_target_gap:
                start_pos = pos
                # Find insertion end position
                while pos < seq_len:
                    # Check if target species still has base
                    target_has_base = species_data[target]['seq'][pos] != '-'
                    
                    # Check if non-target species still have gaps
                    non_target_gap = all(
                        species_data[sp]['seq'][pos] == '-' 
                        for sp in common_species
                    )
                    
                    if not (target_has_base and non_target_gap):
                        break
                    pos += 1
                
                end_pos = pos
                insertion_length = end_pos - start_pos
                
                # Check if insertion length meets minimum requirement
                if insertion_length >= min_size:
                    # Calculate genomic coordinates for target species
                    target_info = species_data[target]
                    target_start, target_end = calc_genomic_interval(
                        target_info['start'],
                        target_info['size'],
                        target_info['strand'],
                        target_info['seq'],
                        start_pos,
                        end_pos,
                        target_info['chrom_size']
                    )
                    
                    # Calculate reference species coordinates (insertion point)
                    ref_info = species_data[ref_species]
                    # Reference species should have gaps in insertion region, so coordinate is a point
                    ref_before = 0
                    for i in range(start_pos):
                        if ref_info['seq'][i] != '-':
                            ref_before += 1
                    
                    if ref_info['strand'] == '+':
                        ref_point = ref_info['start'] + ref_before
                    else:
                        ref_point = ref_info['chrom_size'] - (ref_info['start'] + ref_before)
                    
                    sv_events.append({
                        'type': 'INS',
                        'size': insertion_length,
                        'ref_chrom': ref_info['chrom'],
                        'ref_start': ref_point,
                        'ref_end': ref_point + 1,  # Point position
                        'ref_strand': ref_info['strand'],
                        'target_chrom': target_info['chrom'],
                        'target_start': target_start,
                        'target_end': target_end,
                        'target_strand': target_info['strand'],
                        'block_idx': block_idx,
                        'sv_source': 'intra-block',
                        'gap_type': 'target_insertion'
                    })
            else:
                pos += 1
        
        # ================= 2. Detect target species-specific deletions =================
        # Target species has gap, while most other species have bases
        pos = 0
        while pos < seq_len:
            # Target species must have gap
            target_is_gap = species_data[target]['seq'][pos] == '-'
            
            if target_is_gap:
                start_pos = pos
                # Extend deletion region until condition not met
                while pos < seq_len:
                    target_is_gap = species_data[target]['seq'][pos] == '-'
                    if not target_is_gap:
                        break
                    pos += 1
                end_pos = pos
                deletion_length = end_pos - start_pos
                
                if deletion_length >= min_size:
                    # Check gap ratio of non-target species in this region
                    is_valid = True
                    for sp in common_species:
                        sp_region = species_data[sp]['seq'][start_pos:end_pos]
                        gap_count = sp_region.count('-')
                        gap_ratio = gap_count / deletion_length
                        if gap_ratio > max_gap_ratio:
                            is_valid = False
                            break
                    
                    if not is_valid:
                        continue  # Skip regions that don't meet criteria
                    
                    # Calculate genomic coordinates for target species (deletion point)
                    target_info = species_data[target]
                    # Target species has gaps in deletion region, so coordinate is a point
                    target_before = 0
                    for i in range(start_pos):
                        if target_info['seq'][i] != '-':
                            target_before += 1
                    
                    if target_info['strand'] == '+':
                        target_point = target_info['start'] + target_before
                    else:
                        target_point = target_info['chrom_size'] - (target_info['start'] + target_before)
                    
                    # Calculate reference species coordinates in deletion region
                    ref_info = species_data[ref_species]
                    ref_start, ref_end = calc_genomic_interval(
                        ref_info['start'],
                        ref_info['size'],
                        ref_info['strand'],
                        ref_info['seq'],
                        start_pos,
                        end_pos,
                        ref_info['chrom_size']
                    )
                    
                    sv_events.append({
                        'type': 'DEL',
                        'size': deletion_length,
                        'ref_chrom': ref_info['chrom'],
                        'ref_start': ref_start,
                        'ref_end': ref_end,
                        'ref_strand': ref_info['strand'],
                        'target_chrom': target_info['chrom'],
                        'target_start': target_point,
                        'target_end': target_point + 1,  # Point position
                        'target_strand': target_info['strand'],
                        'block_idx': block_idx,
                        'sv_source': 'intra-block',
                        'gap_type': 'target_deletion'
                    })
            else:
                pos += 1
    
    return sv_events

def detect_inter_block_deletions(blocks, target, ref_species, non_target_species, max_target_gap=10, min_ref_gap=50):
    """Inter-block deletion detection"""
    deletions = []
    target_list = [target]
    
    for i in range(len(blocks)-1):
        prev_block = blocks[i]['species']
        curr_block = blocks[i+1]['species']
        
        # 1. Check if all target species are strictly continuous
        target_valid = True
        for sp in target_list:
            if sp not in prev_block or sp not in curr_block:
                target_valid = False
                break
                
            prev = prev_block[sp]
            curr = curr_block[sp]
            
            # Check chromosome and strand direction consistency
            if prev['chrom'] != curr['chrom'] or prev['strand'] != curr['strand']:
                target_valid = False
                break
                
            gap = calculate_gap(prev, curr, sp)
            if gap is None or abs(gap) > max_target_gap:
                target_valid = False
                break
                
        if not target_valid:
            continue
            
        # 2. Get common species (excluding target species)
        common_species = set(prev_block.keys()) & set(curr_block.keys()) - set(target_list)
        if not common_species:
            continue
            
        # 3. Check if all non-target species must have gaps
        non_target_gaps = {}
        all_non_target_valid = True
        for sp in common_species:
            prev = prev_block[sp]
            curr = curr_block[sp]
            
            # Check chromosome and strand direction consistency
            if prev['chrom'] != curr['chrom'] or prev['strand'] != curr['strand']:
                all_non_target_valid = False
                break
                
            gap = calculate_gap(prev, curr, sp)
            if gap is None or abs(gap) < min_ref_gap:  # Non-target species must be discontinuous
                all_non_target_valid = False
                break
                
            non_target_gaps[sp] = gap
            
        if not all_non_target_valid:
            continue
            
        # Use reference species for coordinates
        if ref_species not in prev_block or ref_species not in curr_block:
            continue
            
        ref_prev = prev_block[ref_species]
        ref_curr = curr_block[ref_species]
        ref_gap = calculate_gap(ref_prev, ref_curr, ref_species)
        if ref_gap is None:
            continue
            
        del_size = abs(ref_gap)
        
        # 5. Coordinate calculation
        ref_start = ref_prev['end']
        ref_end = ref_curr['start']
        
        # Target species position should be point coordinate
        #target_pos = prev_block[target]['actual_end']
        if prev_block[target]['strand'] == '+':
            target_pos = prev_block[target]['actual_end']
        # For negative strand, use chromosome total length minus start and length of previous block
        else:
            target_pos = prev_block[target]['chrom_size'] - (prev_block[target]['start'] + prev_block[target]['size'])
        
        
        deletions.append({
            'type': 'DEL',
            'size': del_size,
            'ref_chrom': ref_prev['chrom'],
            'ref_start': min(ref_start, ref_end),
            'ref_end': max(ref_start, ref_end),
            'ref_strand': ref_prev['strand'],
            'target_chrom': prev_block[target]['chrom'],
            'target_start': target_pos,
            'target_end': target_pos + 1,  # Point position represented as (start, start+1)
            'target_strand': prev_block[target]['strand'],
            'block1_idx': i,
            'block2_idx': i+1,
            'non_target_gaps': non_target_gaps,
            'sv_source': 'inter-block'
        })
    
    return deletions

def detect_inter_block_insertions(blocks, target, ref_species, non_target_species, min_target_gap=50, max_ref_gap=10, translocation_blocks=set()):
    """Distinguish insertions from translocations"""
    insertions = []
    
    for i in range(len(blocks)-1):
        # Skip blocks involved in translocations
        if i in translocation_blocks or (i+1) in translocation_blocks:
            continue
            
        prev_block = blocks[i]['species']
        curr_block = blocks[i+1]['species']
        
        # 1. Target species must have gap
        if target not in prev_block or target not in curr_block:
            continue
            
        t_prev = prev_block[target]
        t_curr = curr_block[target]
        
        # Check chromosome and strand direction consistency
        if t_prev['chrom'] != t_curr['chrom'] or t_prev['strand'] != t_curr['strand']:
            continue
            
        target_gap = calculate_gap(t_prev, t_curr, target)
        if target_gap is None or abs(target_gap) < min_target_gap:
            continue
            
        # 2. Non-target species must be continuous
        valid = True
        for species in non_target_species + [ref_species]:
            if species == target:
                continue
                
            if species not in prev_block or species not in curr_block:
                valid = False
                break
                
            prev = prev_block[species]
            curr = curr_block[species]
            
            # Check chromosome and strand direction consistency
            if prev['chrom'] != curr['chrom'] or prev['strand'] != curr['strand']:
                valid = False
                break
                
            gap = calculate_gap(prev, curr, species)
            if gap is None or abs(gap) > max_ref_gap:
                valid = False
                break
                
        if not valid:
            continue
            
        # 3. Ensure correct target position order
        if t_prev['strand'] == '+':
            target_start = t_prev['actual_end']
            target_end = t_curr['actual_start']
        else:
            target_start = t_curr['actual_end']
            target_end = t_prev['actual_start']
        
        # If order is reversed, it might be translocation, skip
        if target_start >= target_end:
            continue
            
        # 4. Calculate insertion size and check if ≤150000
        target_gap_size = abs(target_gap)
        if target_gap_size > 150000:  # New condition: only accept insertions ≤150000
            continue
            
        # 5. Create insertion record
        ref_pos = prev_block[ref_species]['end']
        
        insertions.append({
            'type': 'INS',
            'size': target_gap_size,
            'ref_chrom': prev_block[ref_species]['chrom'],
            'ref_start': ref_pos,
            'ref_end': ref_pos + 1,  # Point position
            'ref_strand': prev_block[ref_species]['strand'],
            'target_chrom': t_prev['chrom'],
            'target_start': target_start,
            'target_end': target_end,
            'target_strand': t_prev['strand'],
            'block1_idx': i,
            'block2_idx': i+ 1,
            'sv_source': 'inter-block'
        })
    
    return insertions

def detect_translocations(blocks, target, ref_species, non_target_species, max_gap=100):
    """Block1 and block3 chromosomes are same, determine translocation type based on block2 chromosome"""
    translocations = []
    translocation_block_indices = set()
    
    for i in range(len(blocks)-2):
        b1 = blocks[i]['species']
        b2 = blocks[i+1]['species']
        b3 = blocks[i+2]['species']
        
        # 1. Target species must exist in all three blocks
        if target not in b1 or target not in b2 or target not in b3:
            continue
            
        t1 = b1[target]
        t2 = b2[target]
        t3 = b3[target]
        
        # 2. Require block1 and block3 chromosomes to be same, otherwise skip
        if t1['chrom'] != t3['chrom']:
            continue

        # 2.1 All three blocks' target species strand directions must be consistent
        if t1['strand'] != t2['strand'] or t2['strand'] != t3['strand']:
            continue    
            
        # 3. Check if non-target species gaps are within max_gap range
        non_target_valid = True
        for species in non_target_species + [ref_species]:
            if species not in b1 or species not in b2 or species not in b3:
                non_target_valid = False
                break
                
            s1 = b1[species]
            s2 = b2[species]
            s3 = b3[species]
            
            # Check chromosome consistency
            if s1['chrom'] != s2['chrom'] or s2['chrom'] != s3['chrom']:
                non_target_valid = False
                break
                
            # Check strand direction consistency
            if s1['strand'] != s2['strand'] or s2['strand'] != s3['strand']:
                non_target_valid = False
                break
                
            # Check if inter-block gaps are within max_gap range
            gap1 = calculate_gap(s1, s2, species)
            gap2 = calculate_gap(s2, s3, species)
            
            if gap1 is None or gap2 is None or abs(gap1) > max_gap or abs(gap2) > max_gap:
                non_target_valid = False
                break
                
        if not non_target_valid:
            continue
            
        # 4. Check if target species block2 is not between block1 and block3 coordinates
        if t1['strand'] == '+':
            # Positive strand: coordinates increase from small to large
            block1_start = t1['actual_start']
            block2_start = t2['actual_start']
            block2_end = t2['actual_end']
            block3_end = t3['actual_end']
            
            # Check if block2 is not between block1 and block3
            is_translocation = (block2_start < block1_start and block2_end < block1_start) or (block2_start > block3_end and block2_end > block3_end)
            
            # Check gap between block1 and block3
            gap_b1_b3 = t3['actual_start'] - t1['actual_end']
        else:
            # Negative strand: coordinates decrease from large to small
            block1_start = t1['actual_start']
            block2_start = t2['actual_start']
            block2_end = t2['actual_end']
            block3_end = t3['actual_end']
            
            # Check if block2 is not between block1 and block3
            is_translocation = (block2_start < block1_start and block2_end < block1_start) or (block2_start > block3_end and block2_end > block3_end)
            
            # Check gap between block1 and block3
            gap_b1_b3 = t1['actual_end'] - t3['actual_start']
        
        # Ensure gap between target species block1 and block3 is within 100bp range
        if abs(gap_b1_b3) > 100:
            continue
            
        if not is_translocation:
            continue
            
        # 5. Determine translocation type
        if t2['chrom'] == t1['chrom']:
            trans_type = 'INTRACHROM_TRANS'
            original_chrom = t1['chrom']
            original_start = min(t1['actual_start'], t3['actual_start'])
            original_end = max(t1['actual_end'], t3['actual_end'])
        else:
            trans_type = 'INTERCHROM_TRANS'
            original_chrom = t1['chrom']
            original_start = t1['actual_start']
            original_end = t3['actual_end']
        
        # 6. Create translocation record
        ref_info = b1[ref_species] if ref_species in b1 else None
        
        translocations.append({
            'type': trans_type,
            'size': t2['size'],
            'ref_chrom': ref_info['chrom'] if ref_info else 'NA',
            'ref_start': ref_info['end'] if ref_info else 0,
            'ref_end': b3.get(ref_species, {}).get('actual_start', ref_info['end'] if ref_info else 0),
            'ref_strand': ref_info['strand'] if ref_info else 'NA',
            'target_chrom': t2['chrom'],
            'target_start': t2['actual_start'],
            'target_end': t2['actual_end'],
            'target_strand': t2['strand'],
            'target_tra_chrom': t1['chrom'],
            'target_tra_start': t1['actual_end'],
            'target_tra_end': t3['actual_start'] + 1,
            'target_tra_strand': t1['strand'],
            'original_chrom': original_chrom,
            'original_start': original_start,
            'original_end': original_end,
            'block1_idx': i,
            'block2_idx': i+1,
            'block3_idx': i+2,
            'sv_source': 'inter-block',
            'gap_b1_b3': gap_b1_b3  # Record gap between block1 and block3
        })
        translocation_block_indices.add(i+1)
    
    return translocations, translocation_block_indices

def detect_inversions(blocks, target, ref_species, non_target_species, max_gap=100, inv_more_block_max=1000):
    """Detect inversion variants, including three-block inversions and multi-block inversions"""
    inversions = []
    # First detect three-block inversions
    inversions.extend(detect_three_block_inversions(blocks, target, ref_species, non_target_species, max_gap))
    # Detect multi-block inversions (four or more blocks)
    inversions.extend(detect_multi_block_inversions(blocks, target, ref_species, non_target_species, max_gap, inv_more_block_max))
    return inversions

def detect_three_block_inversions(blocks, target, ref_species, non_target_species, max_gap):
    """Detect classic three-block inversion patterns (+-+ or -+-)"""
    three_block_inversions = []
    for i in range(len(blocks)-2):
        b1 = blocks[i]['species']
        b2 = blocks[i+1]['species']
        b3 = blocks[i+2]['species']
        
        # 1. Target species must exist in all three blocks
        if target not in b1 or target not in b2 or target not in b3:
            continue
            
        t1 = b1[target]
        t2 = b2[target]
        t3 = b3[target]
        
        # 2. Check chromosome consistency
        if t1['chrom'] != t2['chrom'] or t2['chrom'] != t3['chrom']:
            continue
            
        # 3. Check strand direction pattern (+-+ or -+-)
        if not ((t1['strand'] == '+' and t2['strand'] == '-' and t3['strand'] == '+') or
                (t1['strand'] == '-' and t2['strand'] == '+' and t3['strand'] == '-')):
            continue
            
        # 4. Check non-target species continuity
        valid = True 
        for species in non_target_species + [ref_species]:
            if species not in b1 or species not in b2 or species not in b3:
                valid = False
                break
                
            s1 = b1[species]
            s2 = b2[species]
            s3 = b3[species]
            
            # Check chromosome consistency
            if s1['chrom'] != s2['chrom'] or s2['chrom'] != s3['chrom']:
                valid = False
                break
                
            # Check strand direction consistency
            if s1['strand'] != s2['strand'] or s2['strand'] != s3['strand']:
                valid = False
                break
                
            # Check inter-block gaps
            gap1 = calculate_gap(s1, s2, species)
            gap2 = calculate_gap(s2, s3, species)
            
            if gap1 is None or gap2 is None or abs(gap1) > max_gap or abs(gap2) > max_gap:
                valid = False
                break
                
        if not valid:
            continue
            
        # 5. Create inversion record
        ref_info = b1[ref_species]
        inversion_length = t2['size']
        
        three_block_inversions.append({
            'type': 'INV',
            'size': inversion_length,
            'ref_chrom': ref_info['chrom'],
            'ref_start': b1[ref_species]['end'],
            'ref_end': b3[ref_species]['start'],
            'ref_strand': ref_info['strand'],
            'target_chrom': t2['chrom'],
            'target_start': t2['actual_start'],
            'target_end': t2['actual_end'],
            'target_strand': t2['strand'],
            'block1_idx': i,
            'block2_idx': i+1,
            'block3_idx': i+2,
            'sv_source': 'inter-block',
            'pattern': f"{t1['strand']}{t2['strand']}{t3['strand']}"
        })
    
    return three_block_inversions

def detect_multi_block_inversions(blocks, target, ref_species, non_target_species, max_gap, inv_more_block_max):
    """Detect multi-block inversion patterns (e.g., +--+, -++- etc.), using inv_more_block_max as gap threshold"""
    multi_block_inversions = []
    
    # Start checking from 4 blocks up to inv_more_block_max blocks
    for block_count in range(4, inv_more_block_max + 1):
        for i in range(len(blocks) - block_count + 1):
            # Get all blocks in current window
            window_blocks = [blocks[i+j] for j in range(block_count)]
            window_species = [block['species'] for block in window_blocks]
            
            # 1. Check if target species exists in all blocks and chromosomes consistent
            target_valid = True
            target_strands = []
            target_chrom = None
            
            for j, block in enumerate(window_species):
                if target not in block:
                    target_valid = False
                    break
                    
                t_data = block[target]
                if target_chrom is None:
                    target_chrom = t_data['chrom']
                elif t_data['chrom'] != target_chrom:
                    target_valid = False
                    break
                    
                target_strands.append(t_data['strand'])
                
            if not target_valid:
                continue
                
            # 2. Check strand direction pattern (e.g., +--+, -++- etc.)
            # First and last strand directions must be same, all middle blocks' strands must be same and opposite to first/last
            expected_middle_strand = '-' if target_strands[0] == '+' else '+'
            pattern_valid = True
            
            # Check first and last
            if target_strands[0] != target_strands[-1]:
                pattern_valid = False
            # Check all middle blocks
            for j in range(1, block_count - 1):
                if target_strands[j] != expected_middle_strand:
                    pattern_valid = False
                    break
                    
            if not pattern_valid:
                continue
                
            # 3. Check non-target species continuity and gaps
            non_target_valid = True
            for species in non_target_species + [ref_species]:
                # Check if species exists in all blocks
                for j, block in enumerate(window_species):
                    if species not in block:
                        non_target_valid = False
                        break
                        
                if not non_target_valid:
                    break
                    
                # Check chromosome consistency
                species_chrom = None
                for j, block in enumerate(window_species):
                    s_data = block[species]
                    if species_chrom is None:
                        species_chrom = s_data['chrom']
                    elif s_data['chrom'] != species_chrom:
                        non_target_valid = False
                        break
                        
                if not non_target_valid:
                    break
                    
                # Check strand direction consistency (all blocks' strands must be same)
                species_strand = window_species[0][species]['strand']
                for j, block in enumerate(window_species):
                    if block[species]['strand'] != species_strand:
                        non_target_valid = False
                        break
                        
                if not non_target_valid:
                    break
                    
                # Check all inter-block gaps - use inv_more_block_max as threshold
                for j in range(len(window_species) - 1):
                    prev = window_species[j][species]
                    curr = window_species[j+1][species]
                    gap = calculate_gap(prev, curr, species)
                    
                    if gap is None or abs(gap) > inv_more_block_max:  # Use inv_more_block_max as gap threshold
                        non_target_valid = False
                        break
                        
                if not non_target_valid:
                    break
                    
            if not non_target_valid:
                continue
                
            # 4. Create multi-block inversion record
            # Inversion region is all middle blocks
            middle_blocks = window_species[1:-1]
            total_inversion_size = sum(block[target]['size'] for block in middle_blocks)
            
            # Use first block's reference species information as coordinate reference
            ref_info = window_species[0][ref_species]
            
            # Target species inversion region start and end positions (considering strand direction)
            if target_strands[0] == '+':
                target_start = window_species[1][target]['actual_start']
                target_end = window_species[-2][target]['actual_end']
            else:
                target_start = window_species[-2][target]['actual_start']
                target_end = window_species[1][target]['actual_end']
                
            # Build pattern string
            pattern = ''.join(target_strands)
            
            multi_block_inversions.append({
                'type': 'INV',
                'size': total_inversion_size,
                'ref_chrom': ref_info['chrom'],
                'ref_start': window_species[0][ref_species]['end'],
                'ref_end': window_species[-1][ref_species]['start'],
                'ref_strand': ref_info['strand'],
                'target_chrom': target_chrom,
                'target_start': target_start,
                'target_end': target_end,
                'target_strand': target_strands[-2],
                'block_indices': list(range(i, i + block_count)),
                'sv_source': 'inter-block',
                'pattern': pattern,
                'block_count': block_count
            })
    
    return multi_block_inversions

def write_bed(output_file, sv_events, sv_type):
    """Write detected SVs to BED file"""
    with open(output_file, 'w') as f:
        header = [
            "#Ref_Chrom", "Ref_Start", "Ref_End", "Ref_Strand",
            "Target_Chrom", "Target_Start", "Target_End", "Target_Strand",
            "Original_tra_Chrom", "Original_tra_Start", "Original_tra_End", "Original_tra_Strand",
            "SV_Type", "SV_Size", "SV_Source",
            "Block_Indexes", "Additional_Info"
        ]
        f.write("\t".join(header) + "\n")
        
        for event in sv_events:
            if 'block_idx' in event:  # Single block SV
                blocks = str(event['block_idx'])
                additional_info = "intra-block" if 'sv_source' in event and event['sv_source'] == 'intra-block' else "isolated-block"
            elif 'block_indices' in event:  # Multi-block inversion event
                blocks = ",".join(str(x) for x in event['block_indices']) # Convert list to comma-separated string
                additional_info = f"multi-block-inversion (pattern: {event.get('pattern', 'N/A')})"
            elif 'block1_idx' in event and 'block2_idx' in event and 'block3_idx' in event:  # Three-block inversion event
                blocks = f"{event['block1_idx']},{event['block2_idx']},{event['block3_idx']}"
                additional_info = "inter-block"
            else:  # Two-block SVs (or other cases)
                blocks = f"{event['block1_idx']},{event['block2_idx']}"
                additional_info = "inter-block"
            
            # Organize fields
            fields = [
                event['ref_chrom'],
                str(event['ref_start']),
                str(event['ref_end']),
                event.get('ref_strand', '+'),
                event['target_chrom'],
                str(event['target_start']),
                str(event['target_end']),
                event.get('target_strand', '+'),
                event.get('target_tra_chrom', 'NA'),       # Use get method with default 'NA'
                str(event.get('target_tra_start', 'NA')),     # Provide default value 0
                str(event.get('target_tra_end', 'NA')),       # Provide default value 0
                event.get('target_tra_strand', 'NA'),       # Provide default value '+'
                event['type'],
                str(event['size']),
                event.get('sv_source', 'unknown'),
                blocks,
                additional_info
            ]
            f.write("\t".join(fields) + "\n")

def main():
    parser = argparse.ArgumentParser(description='Detect structural variants in target species')
    parser.add_argument('--maf', required=True, help='Input MAF file')
    parser.add_argument('--target', required=True, help='Target species name')
    parser.add_argument('--reference', required=True, help='Reference species name')
    parser.add_argument('--non_target', default="", help='Other non-target species (comma-separated)')
    parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    
    # Optional parameters
    parser.add_argument('--max_target_gap', type=int, default=10, help='Maximum continuous gap for target species, used in deletion detection')
    parser.add_argument('--min_ref_gap', type=int, default=50, help='Minimum SV gap for reference species')
    parser.add_argument('--min_target_gap', type=int, default=50, help='Minimum SV gap for target species')
    parser.add_argument('--max_ref_gap', type=int, default=10, help='Maximum continuous gap for reference species, used in insertion detection')
    parser.add_argument('--max_gap', type=int, default=100, help='Maximum gap for continuous regions, used in translocation and inversion detection')
    parser.add_argument('--inv_more_block_max', type=int, default=1000, help='Maximum allowed gap threshold for all species in multi-block inversion detection')
    parser.add_argument('--min_intra_size', type=int, default=10, help='Minimum size for intra-block SVs')
    parser.add_argument('--max_gap_ratio', type=float, default=0.2, help='Maximum allowed gap ratio for non-target species in deletion detection')
    
    args = parser.parse_args()
    
    # Parse non-target species
    non_target_species = args.non_target.split(',') if args.non_target else []
    
    # Read MAF file
    print(f"Parsing MAF file: {args.maf}")
    blocks = parse_maf_blocks(args.maf)
    print(f"Parsing completed, found {len(blocks)} alignment blocks")
    
    # Detect various SVs
    print("\nDetecting structural variants...")
    
    # First detect translocation events and get translocation block indices
    translocations, translocation_block_indices = detect_translocations(
        blocks, args.target, args.reference, non_target_species, args.max_gap
    )
    
    # Separate two translocation types
    intra_trans = [t for t in translocations if t['type'] == 'INTRACHROM_TRANS']
    inter_trans = [t for t in translocations if t['type'] == 'INTERCHROM_TRANS']
    
    # Intra-block SV detection
    intra_block_svs = detect_block_svs(
        blocks, args.target, args.reference, non_target_species,
        args.min_intra_size, args.max_gap_ratio
    )
    
    # Inter-block deletion detection
    inter_block_deletions = detect_inter_block_deletions(
        blocks, args.target, args.reference, non_target_species,
        args.max_target_gap, args.min_ref_gap
    )
    
    # Inter-block insertion detection (pass translocation block indices)
    inter_block_insertions = detect_inter_block_insertions(
        blocks, args.target, args.reference, non_target_species,
        args.min_target_gap, args.max_ref_gap, translocation_block_indices
    )
    
    # Inversion detection (including three-block and multi-block)
    inversions = detect_inversions(
        blocks, args.target, args.reference, non_target_species,
        args.max_gap, args.inv_more_block_max
    )
    
    # Separate intra-block insertions and deletions
    intra_ins = [sv for sv in intra_block_svs if sv['type'] == 'INS']
    intra_del = [sv for sv in intra_block_svs if sv['type'] == 'DEL']
    
    # Merge all SVs
    all_deletions = inter_block_deletions + intra_del
    all_insertions = inter_block_insertions + intra_ins
    
    # Create output directory
    os.makedirs(os.path.dirname(args.output_prefix) or '.', exist_ok=True)
    
    # Write result files
    write_bed(f"{args.output_prefix}_deletions.bed", all_deletions, "DEL")
    write_bed(f"{args.output_prefix}_insertions.bed", all_insertions, "INS")
    write_bed(f"{args.output_prefix}_inversions.bed", inversions, "INV")
    write_bed(f"{args.output_prefix}_intrachrom_translocations.bed", intra_trans, "TRANS")
    write_bed(f"{args.output_prefix}_interchrom_translocations.bed", inter_trans, "TRANS")
    
    # Output statistics
    print("\nDetection results statistics:")
    print(f"Deletions (DEL): {len(all_deletions)} (intra-block: {len(intra_del)}, inter-block: {len(inter_block_deletions)})")
    print(f"Insertions (INS): {len(all_insertions)} (intra-block: {len(intra_ins)}, inter-block: {len(inter_block_insertions)})")
    print(f"Inversions (INV): {len(inversions)}")
    print(f"Intrachromosomal translocations: {len(intra_trans)}")
    print(f"Interchromosomal translocations: {len(inter_trans)}")
    print(f"Total SV events: {len(all_deletions)+len(all_insertions)+len(inversions)+len(translocations)}")
    
    # Output file paths
    print("\nResult files:")
    print(f"Deletions: {args.output_prefix}_deletions.bed")
    print(f"Insertions: {args.output_prefix}_insertions.bed")
    print(f"Inversions: {args.output_prefix}_inversions.bed")
    print(f"Intrachromosomal translocations: {args.output_prefix}_intrachrom_translocations.bed")
    print(f"Interchromosomal translocations: {args.output_prefix}_interchrom_translocations.bed")

if __name__ == "__main__":
    main()