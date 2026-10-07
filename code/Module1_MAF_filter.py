# -*- coding: utf-8 -*- 
import argparse

def parse_maf(maf_file): 
    current_block = [] 
    with open(maf_file, 'r') as f: 
        for line in f: 
            line = line.strip() 
            if not line or line.startswith('#'): 
                continue
            if line.startswith('a'): 
                if current_block: 
                    yield current_block
                current_block = [line] 
            else: 
                current_block.append(line) 
    if current_block: 
        yield current_block

def calculate_gap_ratio(sequences, start, end): 
    total = (end - start) * len(sequences) 
    if total == 0: 
        return 0.0
    gaps = sum(seq[start:end].count('-') for seq in sequences) 
    return gaps / total

def parse_thresholds(s): 
    return {k: float(v) for pair in s.split(',') for k, v in [pair.split(':')]} if s else {} 

def filter_block(block, args): 
    species_data = {} 
    ref_seq, ref_len = None, 0
    required = set(args.required_species.split(',')) if args.required_species else set() 
    species_gap = parse_thresholds(args.species_gap_ratio) 
    
    score_val = None
    for line in block: 
        if line.startswith('a'): 
            parts = line.split() 
            for part in parts: 
                if part.startswith('score='): 
                    try: 
                        score_val = float(part.split('=')[1]) 
                    except ValueError: 
                        pass 
            break

    # Parse sequences
    for line in block: 
        if line.startswith('s'): 
            parts = line.split() 
            species = parts[1].split('.')[0] 
            seq = parts[-1].upper() 
            species_data[species] = seq
            if species == args.ref_species: 
                ref_seq = seq
                ref_len = int(parts[3])  # MAF alignment length
    
    # === Basic Information ===
    print(f"\n[DEBUG] Processing new block | Reference species={args.ref_species} | Score={score_val} | Alignment length={ref_len} | Number of species={len(species_data)}") 
    
    if args.abnormal == 'yes' and score_val is not None:
        if score_val < 0:
            print(f"  → Filter reason: Negative score ({score_val})") 
            return False

    # Condition 1: Required species must be present
    if required and not required.issubset(species_data.keys()): 
        print(f"  → Filter reason: Missing required species {required - species_data.keys()}") 
        return False
    
    # Condition 2: Reference species existence
    if not ref_seq: 
        print("  → Filter reason: Reference species not found") 
        return False
    
    # Condition 3: Minimum number of sequences
    if len(species_data) < args.min_sequences: 
        print(f"  → Filter reason: Insufficient sequences ({len(species_data)} < {args.min_sequences})") 
        return False
    
    # Condition 4: Reference sequence alignment length
    if ref_len < args.min_ref_length: 
        print(f"  → Filter reason: Alignment length insufficient ({ref_len} < {args.min_ref_length})") 
        return False
    
    # === Global Gap Detection ===
    actual_seq_len = len(ref_seq) 
    region = min(args.gap_region_length, actual_seq_len) 
    left_ratio = calculate_gap_ratio(species_data.values(), 0, region) 
    right_ratio = calculate_gap_ratio(species_data.values(), max(0, actual_seq_len - region), actual_seq_len) 
    
    # === Debug Output 2: Global Gaps ===
    print(f"  → Global GAP detection | Left region:0-{region} | Left gap ratio={left_ratio:.2f} | Right region:{max(0, actual_seq_len - region)}-{actual_seq_len} | Right gap ratio={right_ratio:.2f}") 
    
    if left_ratio > args.max_gap_ratio or right_ratio > args.max_gap_ratio: 
        print(f"  → Filter reason: Global gap exceeds threshold (threshold={args.max_gap_ratio})") 
        return False
    
    # === Species-specific Detection ===
    if species_gap: 
        check_len = args.species_gap_region_length or args.gap_region_length
        for species, max_ratio in species_gap.items(): 
            if species not in species_data: 
                continue
            seq = species_data[species] 
            slen = len(seq) 
            if slen == 0: 
                continue
            
            # Left end detection
            left_end = min(check_len, slen) 
            left = seq[:left_end].count('-') / left_end
            
            # Right end detection
            right_start = max(0, slen - check_len) 
            right = seq[right_start:].count('-') / (slen - right_start) if (slen - right_start) >0 else 0
            
            # === Debug Output 3: Species Detection ===
            print(f"  → Species detection [{species}] | Left gap={left:.2f} (threshold={max_ratio}) | Right gap={right:.2f} (threshold={max_ratio})") 
            
            if left > max_ratio or right > max_ratio: 
                print(f"  → Filter reason: Species {species} gap exceeds threshold") 
                return False
    
    print("  → Block retained") 
    return True

def main(): 
    parser = argparse.ArgumentParser(description='MAF filtering tool') 
    parser.add_argument('--input', required=True) 
    parser.add_argument('--output', required=True) 
    parser.add_argument('--min_sequences', type=int, default=2) 
    parser.add_argument('--required_species', type=str) 
    parser.add_argument('--ref_species', default='ho_sa') 
    parser.add_argument('--min_ref_length', type=int, default=100) 
    parser.add_argument('--max_gap_ratio', type=float, default=0.05) 
    parser.add_argument('--gap_region_length', type=int, default=100) 
    parser.add_argument('--species_gap_ratio', type=str) 
    parser.add_argument('--species_gap_region_length', type=int) 
    parser.add_argument('--abnormal', choices=['yes', 'no'], default='yes', 
                        help='Whether to remove blocks with negative scores (yes/no). Default is no.') 
    
    args = parser.parse_args() 

    with open(args.output, 'w') as f: 
        for i, block in enumerate(parse_maf(args.input)): 
            print(f"\n=== Processing block {i+1} ===") 
            if filter_block(block, args): 
                f.write('\n'.join(block) + '\n\n') 

if __name__ == '__main__': 
    main()