import argparse
import os
import sys
from collections import defaultdict

def read_fasta(file_path):
    """Read FASTA file and return sequence dictionary {sequence name: sequence string}"""
    sequences = {}
    with open(file_path, 'r') as f:
        current_id = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    return sequences

def write_fasta(file_path, sequences):
    """Write sequence dictionary to FASTA file"""
    with open(file_path, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

def parse_bed_file(bed_path):
    """Parse enhanced BED file, support any number of insertion commands"""
    events = []
    with open(bed_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):  # Skip empty lines and comment lines
                continue
                
            parts = line.split('\t')
            if len(parts) < 3:
                print(f"Warning: Insufficient columns ({len(parts)} columns) in line {line_num}, skipping: {line}")
                continue
                
            try:
                filename = parts[0]
                num_seqs = int(parts[1])
                
                regions = []
                for i in range(2, 2 + num_seqs * 2, 2):
                    if i >= len(parts):
                        break
                    start = int(parts[i]) if parts[i] != 'NA' else None
                    end = int(parts[i+1]) if i+1 < len(parts) and parts[i+1] != 'NA' else None
                    regions.append((start, end))
                
                insert_targets = []
                for i in range(2 + num_seqs * 2, len(parts)):
                    target_spec = parts[i]
                    if target_spec == 'NA':
                        continue
                    target_parts = target_spec.split(':')
                    if len(target_parts) != 4:
                        print(f"Warning: Invalid insertion target format in line {line_num}: {target_spec}")
                        continue
                    seq_name = target_parts[0]
                    target_file = target_parts[1]
                    chrom = target_parts[2]
                    pos = int(target_parts[3])
                    insert_targets.append((seq_name, target_file, chrom, pos))
                
                if len(regions) != num_seqs:
                    print(f"Warning: Number of region coordinates ({len(regions)}) does not match sequence count ({num_seqs}) in line {line_num}")
                
                event = {
                    'filename': filename,
                    'num_seqs': num_seqs,
                    'regions': regions,
                    'insert_targets': insert_targets,
                    'line_number': line_num
                }
                
                events.append(event)
                print(f"Parsed event {line_num}: {filename} - {num_seqs} sequences, {len(insert_targets)} insertion targets")
                
            except Exception as e:
                print(f"Error parsing BED line {line_num}: {str(e)} - Content: {line}")
    
    return events

def process_insertions(args):
    """Process all insertion operations - Insertions based on original positions, coordinate records consider offsets"""
    print("\nParsing BED file...")
    events = parse_bed_file(args.bed_file)
    total_events = len(events)
    print(f"Found {total_events} insertion events")
    
    seq_data = {}
    target_files = set()
    
    # Step 1: Load all alignment files
    for event in events:
        filename = event['filename']
        if filename in seq_data:
            continue
            
        aln_path = os.path.join(args.aln_dir, filename)
        if not os.path.exists(aln_path):
            print(f"Error: File {aln_path} does not exist")
            continue
        
        sequences = read_fasta(aln_path)
        if not sequences:
            print(f"Error: No sequences in alignment file {aln_path}")
            continue
        
        seq_data[filename] = {
            'sequences': sequences,
            'seq_order': list(sequences.keys())
        }
        print(f"Loaded alignment file {filename}: {len(sequences)} sequences")
    
    # Step 2: Build insertion event list
    all_insert_events = []
    for event in events:
        filename = event['filename']
        if filename not in seq_data:
            continue
        
        sequences = seq_data[filename]['sequences']
        seq_order = seq_data[filename]['seq_order']
        
        for insert_target in event['insert_targets']:
            seq_name, target_file, chrom, pos = insert_target
            
            if seq_name not in sequences:
                matched = False
                for seq_id in sequences.keys():
                    if seq_name in seq_id or seq_id in seq_name:
                        seq_name = seq_id
                        matched = True
                        break
                if not matched:
                    print(f"Warning: Sequence name '{seq_name}' in event {event['line_number']} does not exist in file {filename}, skipping")
                    continue
            
            seq_content = sequences[seq_name]
            seq_len = len(seq_content)
            
            region_start, region_end = None, None
            try:
                seq_index = seq_order.index(seq_name)
                if seq_index < len(event['regions']):
                    region_start, region_end = event['regions'][seq_index]
            except ValueError:
                pass
            
            if region_start is None or region_end is None:
                if event['regions']:
                    region_start, region_end = event['regions'][0]
                    print(f"Warning: Using first region coordinates as default for sequence '{seq_name}'")
            
            ins_event = {
                'source_file': filename,
                'seq_name': seq_name,
                'seq_content': seq_content,
                'seq_len': seq_len,
                'region': (region_start, region_end),
                'target_file': target_file,
                'chrom': chrom,
                'insert_pos': pos,  # Absolute position in original sequence
                'line_number': event['line_number']
            }
            all_insert_events.append(ins_event)
            target_files.add(target_file)
    
    print(f"Generated {len(all_insert_events)} insertion points in total")
    
    # Step 3: Load all target files
    print("\nLoading target sequence files...")
    target_sequences = {}
    for target_file in target_files:
        file_path = os.path.join(args.target_dir, target_file)
        if not os.path.exists(file_path):
            print(f"Error: Target file {file_path} does not exist")
            continue
        target_sequences[target_file] = read_fasta(file_path)
        print(f"Loaded target file {target_file}: {len(target_sequences[target_file])} chromosomes")
    
    # Step 4: Create copies for each target file
    updated_sequences = {}
    for target_file, chrom_seqs in target_sequences.items():
        updated_sequences[target_file] = {chrom: seq for chrom, seq in chrom_seqs.items()}
    
    # Step 5: Organize insertion events: Group by (target_file, chrom)
    insertion_groups = defaultdict(list)
    for event in all_insert_events:
        key = (event['target_file'], event['chrom'])
        insertion_groups[key].append(event)
    
    # Step 6: Process insertion events and record coordinates
    coord_records = []
    for (target_file, chrom), events in insertion_groups.items():
        if target_file not in updated_sequences or chrom not in updated_sequences[target_file]:
            print(f"Warning: Chromosome {chrom} in target file {target_file} does not exist, skipping this group of events")
            continue
        
        # Get original sequence
        original_seq = updated_sequences[target_file][chrom]
        current_seq = original_seq
        
        # Sort by insertion position (calculate offsets from smallest to largest)
        events_sorted_by_pos = sorted(events, key=lambda x: x['insert_pos'])
        
        # Calculate offset for each insertion point (sum of lengths of all previous inserted sequences)
        offset_map = {}
        cumulative_offset = 0
        last_pos = -1
        
        for event in events_sorted_by_pos:
            # Special handling required if multiple insertions at the same position
            if event['insert_pos'] == last_pos:
                cumulative_offset += event['seq_len']
            else:
                cumulative_offset = sum(e['seq_len'] for e in events_sorted_by_pos 
                                       if e['insert_pos'] < event['insert_pos'])
            
            offset_map[(event['source_file'], event['seq_name'])] = cumulative_offset
            last_pos = event['insert_pos']
        
        # Sort events by insertion position from largest to smallest (insert from back to front)
        sorted_events = sorted(events, key=lambda x: x['insert_pos'], reverse=True)
        
        print(f"\nProcessing chromosome {chrom} of target file {target_file}, total {len(events)} insertion events")
        print(f"Original length: {len(original_seq)}")
        
        # Perform insertion operations (based on original positions)
        for event in sorted_events:
            insert_pos = event['insert_pos']
            seq_content = event['seq_content']
            seq_len = event['seq_len']
            
            # Perform insertion (using original position)
            current_seq = current_seq[:insert_pos] + seq_content + current_seq[insert_pos:]
            
            # Calculate new coordinates of special regions (considering lengths of previous inserted sequences)
            offset = offset_map.get((event['source_file'], event['seq_name']), 0)
            new_region_start = None
            new_region_end = None
            if event['region'][0] is not None and event['region'][1] is not None:
                new_region_start = insert_pos + offset + event['region'][0]
                new_region_end = insert_pos + offset + event['region'][1]
            
            # Record coordinates
            record = {
                'source_file': event['source_file'],
                'seq_name': event['seq_name'],
                'target_file': target_file,
                'chrom': chrom,
                'region_start': event['region'][0],
                'region_end': event['region'][1],
                'new_region_start': new_region_start,
                'new_region_end': new_region_end,
                'original_insert_pos': insert_pos,
                'actual_insert_pos': insert_pos,  # Keep insertion position as original
                'offset': offset,  # Offset for coordinate recording
                'seq_len': seq_len
            }
            coord_records.append(record)
            
            print(f"  Insertion event (line {event['line_number']}): Sequence {event['seq_name']} inserted at position {insert_pos}")
            print(f"    Insertion length: {seq_len}, current chromosome length: {len(current_seq)}")
            if new_region_start is not None:
                print(f"    New coordinates of special region: {new_region_start}-{new_region_end} (Original region: {event['region'][0]}-{event['region'][1]}, Offset: {offset})")
        
        # Update sequence
        updated_sequences[target_file][chrom] = current_seq
        print(f"Final length: {len(current_seq)}")
    
    # Step 7: Save results
    print("\nSaving result files...")
    os.makedirs(args.output_dir, exist_ok=True)
    
    for target_file, chroms in updated_sequences.items():
        output_path = os.path.join(args.output_dir, f"updated_{target_file}")
        write_fasta(output_path, chroms)
        print(f"Saved updated target file: {output_path}")
    
    # Modify here: Correctly handle the path of coord_file
    if os.path.isabs(args.coord_file):
        # Use directly if absolute path is provided
        coord_file_path = args.coord_file
    else:
        # Concatenate with output directory if relative path is provided
        coord_file_path = os.path.join(args.output_dir, args.coord_file)
    
    # Ensure the directory of the coordinate file exists
    coord_dir = os.path.dirname(coord_file_path)
    if coord_dir and not os.path.exists(coord_dir):
        os.makedirs(coord_dir, exist_ok=True)
    
    with open(coord_file_path, 'w') as f:
        f.write("SourceFile\tSequenceName\tTargetFile\tChromosome\tRegionStart\tRegionEnd\tNewRegionStart\tNewRegionEnd\tOriginalInsertPos\tActualInsertPos\tOffset\tSequenceLength\n")
        for record in coord_records:
            f.write(f"{record['source_file']}\t")
            f.write(f"{record['seq_name']}\t")
            f.write(f"{record['target_file']}\t")
            f.write(f"{record['chrom']}\t")
            f.write(f"{record['region_start'] or 'NA'}\t")
            f.write(f"{record['region_end'] or 'NA'}\t")
            f.write(f"{record['new_region_start'] or 'NA'}\t")
            f.write(f"{record['new_region_end'] or 'NA'}\t")
            f.write(f"{record['original_insert_pos']}\t")
            f.write(f"{record['actual_insert_pos']}\t")
            f.write(f"{record['offset']}\t")  # New offset field
            f.write(f"{record['seq_len']}\n")
    
    print(f"Coordinate records saved to: {coord_file_path}")
    print("Processing completed!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='DNA Sequence Insertion Tool (Insertions based on original positions, coordinate records consider offsets)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--aln-dir', required=True, 
                        help='Directory containing alignment sequence files')
    parser.add_argument('--target-dir', required=True, 
                        help='Directory containing target genome files')
    parser.add_argument('--bed-file', required=True, 
                        help='BED format coordinate file (no header)')
    parser.add_argument('--output-dir', required=True, 
                        help='Output directory')
    parser.add_argument('--coord-file', required=True, 
                        help='Coordinate record output filename')
    
    args = parser.parse_args()
    
    for dir_path in [args.aln_dir, args.target_dir]:
        if not os.path.exists(dir_path):
            print(f"Error: Directory {dir_path} does not exist")
            sys.exit(1)
    if not os.path.exists(args.bed_file):
        print(f"Error: BED file {args.bed_file} does not exist")
        sys.exit(1)
    
    print("Starting to process all insertion operations...")
    process_insertions(args)