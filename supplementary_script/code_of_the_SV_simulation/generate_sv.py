#!/usr/bin/env python3
import os
import glob
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def biological_inversion(sequence, start, end):
    """Perform biologically accurate inversion: reverse complement the specified region"""
    region = sequence[start-1:end]
    region_seq = Seq(region)
    inverted_region = region_seq.reverse_complement()
    return sequence[:start-1] + str(inverted_region) + sequence[end:]

def extract_and_remove(sequence, regions):
    """Extract specified regions from the sequence and remove these regions"""
    # Sort regions by start position from smallest to largest
    regions_sorted = sorted(regions, key=lambda x: x[0])
    
    # Extract fragments (in original order)
    fragments = [sequence[start-1:end] for start, end, _ in regions_sorted]
    
    # Sort regions by start position from largest to smallest (for deletion)
    regions_sorted_reverse = sorted(regions, key=lambda x: x[0], reverse=True)
    
    # Remove regions
    modified_seq = sequence
    for start, end, _ in regions_sorted_reverse:
        modified_seq = modified_seq[:start-1] + modified_seq[end:]
    
    return modified_seq, ''.join(fragments)

def process_operation(seq_record, operation_type, regions):
    """Process the sequence according to the operation type"""
    seq_str = str(seq_record.seq)
    original_length = len(seq_str)
    extracted_sequence = None
    log_messages = []
    
    # Sort regions by position from largest to smallest (avoid position offset)
    sorted_regions = sorted(regions, key=lambda x: x[0], reverse=True)
    
    if operation_type == "deletion":
        for i, (start, end, _) in enumerate(sorted_regions, 1):
            if start < 1 or end > original_length:
                log_messages.append(f"Error: Deletion region {i} ({start}-{end}) is out of sequence range")
                continue
                
            seq_str = seq_str[:start-1] + seq_str[end:]
            log_messages.append(f"Deleted region {i}: {start}-{end} | New length: {len(seq_str)}")
        
        seq_record.seq = Seq(seq_str)
    
    elif operation_type == "inversion":
        for i, (start, end, _) in enumerate(sorted_regions, 1):
            if start < 1 or end > original_length:
                log_messages.append(f"Error: Inversion region {i} ({start}-{end}) is out of sequence range")
                continue
                
            original_region = seq_str[start-1:end]
            seq_str = biological_inversion(seq_str, start, end)
            inverted_region = seq_str[start-1:start-1+(end-start+1)]
            log_messages.append(f"Inverted region {i}: {start}-{end} | Original: {original_region[:10]}... | After inversion: {inverted_region[:10]}...")
        
        seq_record.seq = Seq(seq_str)
    
    elif operation_type == "extraction":
        modified_seq, extracted_sequence = extract_and_remove(seq_str, regions)
        seq_record.seq = Seq(modified_seq)
        
        for i, (start, end, _) in enumerate(regions, 1):
            region_length = end - start + 1
            log_messages.append(f"Extracted region {i}: {start}-{end} ({region_length} bp)")
        
        log_messages.append(f"Total extracted length: {len(extracted_sequence)} bp")
    
    else:
        log_messages.append(f"Unknown operation type: {operation_type}")
    
    return seq_record, extracted_sequence, log_messages

def parse_bed_file(bed_file):
    """Parse enhanced BED file"""
    operations = {}
    print(f"Parsing BED file: {bed_file}")
    
    with open(bed_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            parts = line.strip().split()
            if len(parts) < 5:
                print(f"Warning: Skipping line {line_num} - Insufficient columns: {line.strip()}")
                continue
            
            filename = parts[0]
            seq_id = parts[1]
            operation_type = parts[4].lower()
            
            try:
                start = int(parts[2])
                end = int(parts[3])
                length = end - start + 1
            except ValueError:
                print(f"Error: Line {line_num} contains non-numeric values: {line.strip()}")
                continue
            
            # Create file entry if it does not exist
            if filename not in operations:
                operations[filename] = {}
            
            # Create sequence entry if it does not exist
            if seq_id not in operations[filename]:
                operations[filename][seq_id] = {"regions": [], "operation": operation_type}
            
            # Add region
            operations[filename][seq_id]["regions"].append((start, end, length))
    
    print(f"Successfully parsed operation information for {len(operations)} files")
    return operations

def process_files_with_bed(bed_file, input_dir, output_dir="processed"):
    """Process all FASTA files according to the BED file"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse BED file
    operations = parse_bed_file(bed_file)
    
    processed_count = 0
    skipped_files = []
    
    # Process each file
    for filename, seq_operations in operations.items():
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)
        
        print(f"\nProcessing file: {filename}")
        
        # Check if the file exists
        if not os.path.exists(input_path):
            print(f"Error: File does not exist - {input_path}")
            skipped_files.append(filename)
            continue
        
        # Read FASTA file
        try:
            records = list(SeqIO.parse(input_path, "fasta"))
        except Exception as e:
            print(f"Error parsing {filename}: {str(e)}")
            skipped_files.append(filename)
            continue
        
        if not records:
            print(f"No sequences in file {filename}, skipping processing")
            skipped_files.append(filename)
            continue
        
        # Display sequence IDs in the file
        seq_ids = [rec.id for rec in records]
        print(f"File contains sequences: {', '.join(seq_ids)}")
        
        # Process operations for each sequence
        modified = False
        new_records = []
        extraction_info = {}
        
        # First process all existing sequences
        for record in records:
            seq_id = record.id
            
            if seq_id in seq_operations:
                op_data = seq_operations[seq_id]
                operation_type = op_data["operation"]
                regions = op_data["regions"]
                
                print(f"  > Processing sequence '{seq_id}' - Operation: {operation_type.upper()}")
                
                # Execute operation
                modified_record, extracted_seq, log_msgs = process_operation(
                    record, operation_type, regions
                )
                
                # Output log
                for msg in log_msgs:
                    print(f"    {msg}")
                
                # Update record
                new_records.append(modified_record)
                modified = True
                
                # Save extracted sequence
                if operation_type == "extraction" and extracted_seq:
                    extraction_info[seq_id] = extracted_seq
            else:
                # No operation, keep original sequence
                new_records.append(record)
        
        # Add extracted sequences as new records
        for seq_id, extracted_seq in extraction_info.items():
            new_record = SeqRecord(
                Seq(extracted_seq),
                id=f"{seq_id}_fragments",
                description=f"extracted_regions={len(extraction_info[seq_id])}"
            )
            new_records.append(new_record)
            print(f"  Added extracted sequence: {new_record.id} ({len(extracted_seq)} bp)")
        
        # Write to output file
        SeqIO.write(new_records, output_path, "fasta")
        print(f"  Processing completed → {output_path}")
        processed_count += 1
    
    # Output summary
    print("\n" + "="*50)
    print("Processing completed!")
    print(f"Successfully processed: {processed_count} files")
    
    if skipped_files:
        print(f"Skipped files: {len(skipped_files)}")
        print("\n".join(skipped_files))
    
    print(f"Output directory: {output_dir}")
    print("="*50)

if __name__ == "__main__":
    # Set command line parameters
    parser = argparse.ArgumentParser(
        description='Integrated Genome Editing Tool - Supports deletion, inversion and extraction operations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-b', '--bed', required=True, 
                       help='Path to BED file (Format: filename sequenceID start end operation_type)')
    parser.add_argument('-i', '--input', required=True, 
                       help='Input directory path (contains FASTA files)')
    parser.add_argument('-o', '--output', default='processed', 
                       help='Output directory (default: processed/)')
    
    args = parser.parse_args()
    
    # Print title
    print("="*70)
    print("Integrated Genome Editing Tool")
    print("Supported operation types: DELETION, INVERSION, EXTRACTION")
    print("="*70)
    print(f"BED file: {args.bed}")
    print(f"Input directory: {args.input}")
    print(f"Output directory: {args.output}")
    print("-"*70)
    
    # Execute processing
    process_files_with_bed(
        bed_file=args.bed,
        input_dir=args.input,
        output_dir=args.output
    )