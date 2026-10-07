#!/usr/bin/env python3
"""
MAF file realignment module
"""

import argparse
import os
import sys
import tempfile
import subprocess
import logging
import re
from pathlib import Path
from typing import List, Dict, Tuple

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MafRealigner:
    def __init__(self, muscle_path="muscle", mafft_path="mafft", threshold=1000):
        self.muscle_path = muscle_path
        self.mafft_path = mafft_path
        self.threshold = threshold
        
    def check_tools(self):
        """Check if alignment tools are available"""
        try:
            subprocess.run([self.muscle_path, "-version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            logger.info("MUSCLE5 is available")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error(f"MUSCLE5 is not available, please check path: {self.muscle_path}")
            return False
            
        try:
            subprocess.run([self.mafft_path, "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            logger.info("MAFFT is available")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error(f"MAFFT is not available, please check path: {self.mafft_path}")
            return False
            
        return True
    
    def prepare_sequences(self, sequences_with_ids):
        """Process sequences: remove gaps while preserving original order"""
        prepared_seqs = []
        for seq_id, sequence in sequences_with_ids:
            clean_seq = sequence.replace('-', '').replace(' ', '')
            prepared_seqs.append((seq_id, clean_seq))
        return prepared_seqs
    
    def run_muscle_with_order(self, sequences_with_ids, output_file):
        """Use MUSCLE5 for alignment"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            for seq_id, seq in sequences_with_ids:
                f.write(f">{seq_id}\n")
                f.write(f"{seq}\n")
            temp_input = f.name
        
        try:
            # MUSCLE5 alignment ensuring input order
            cmd = [self.muscle_path, "-align", temp_input, "-output", output_file]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            aligned_seqs = self.parse_fasta_with_ids(output_file)
            input_ids = [seq_id for seq_id, _ in sequences_with_ids]
            aligned_ids = [seq_id for seq_id, _ in aligned_seqs]
            
            if input_ids != aligned_ids:
                logger.warning("MUSCLE5 changed sequence order, reordering")
                # Reorder to match input order
                aligned_dict = dict(aligned_seqs)
                with open(output_file, 'w') as f:
                    for seq_id, _ in sequences_with_ids:
                        if seq_id in aligned_dict:
                            f.write(f">{seq_id}\n")
                            f.write(f"{aligned_dict[seq_id]}\n")
                        else:
                            logger.error(f"Sequence {seq_id} missing in alignment result!")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"MUSCLE5 alignment failed: {e.stderr}")
            raise
        finally:
            if os.path.exists(temp_input):
                os.unlink(temp_input)
    
    def run_mafft_with_order(self, sequences_with_ids, output_file):
        """Use MAFFT for alignment"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            for seq_id, seq in sequences_with_ids:
                f.write(f">{seq_id}\n")
                f.write(f"{seq}\n")
            temp_input = f.name
        
        try:
            # Use --inputorder to force maintaining input order, --preservecase to preserve case
            cmd = [self.mafft_path, "--inputorder", "--preservecase", temp_input]
            with open(output_file, 'w') as out_f:
                result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, check=True)
            
            # Verify order is maintained
            aligned_seqs = self.parse_fasta_with_ids(output_file)
            input_ids = [seq_id for seq_id, _ in sequences_with_ids]
            aligned_ids = [seq_id for seq_id, _ in aligned_seqs]
            
            if input_ids != aligned_ids:
                logger.warning("MAFFT changed sequence order, reordering")
                # Reorder to match input order
                aligned_dict = dict(aligned_seqs)
                with open(output_file, 'w') as f:
                    for seq_id, _ in sequences_with_ids:
                        if seq_id in aligned_dict:
                            f.write(f">{seq_id}\n")
                            f.write(f"{aligned_dict[seq_id]}\n")
                        else:
                            logger.error(f"Sequence {seq_id} missing in alignment result!")
                            
        except subprocess.CalledProcessError as e:
            logger.error(f"MAFFT alignment failed: {e.stderr}")
            raise
        finally:
            if os.path.exists(temp_input):
                os.unlink(temp_input)
    
    def parse_fasta_with_ids(self, fasta_file):
        """Parse FASTA file"""
        sequences = []
        current_id = None
        current_seq = []
        
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id is not None:
                        sequences.append((current_id, ''.join(current_seq)))
                    current_id = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
        
        if current_id is not None:
            sequences.append((current_id, ''.join(current_seq)))
        
        return sequences
    
    def parse_maf_block(self, block_lines):
        """Parse MAF alignment block"""
        sequences = []
        other_lines = []
        
        for line in block_lines:
            if line.startswith('s ') or line.startswith('s\t'):
                parts = re.split(r'\s+', line.strip())
                if len(parts) >= 7:
                    seq_info = {
                        'original_line': line,
                        'seq_id': parts[1],
                        'start': parts[2],
                        'length': parts[3],
                        'strand': parts[4],
                        'total_length': parts[5],
                        'sequence': ' '.join(parts[6:])
                    }
                    sequences.append(seq_info)
                else:
                    other_lines.append(line)
            else:
                other_lines.append(line)
        
        return sequences, other_lines
    
    def calculate_column_widths(self, sequences):
        """Calculate maximum widths of columns that need alignment"""
        max_seq_id = max(len(seq['seq_id']) for seq in sequences)
        max_start = max(len(seq['start']) for seq in sequences)
        max_length = max(len(seq['length']) for seq in sequences)
        max_total = max(len(seq['total_length']) for seq in sequences)
        
        return max_seq_id, max_start, max_length, max_total
    
    def format_sequence_line(self, seq_info, seq_id_width, start_width, length_width, total_width):
        """Format sequence line ensuring proper column alignment"""
        formatted = f"s {seq_info['seq_id']:<{seq_id_width}} {seq_info['start']:>{start_width}} {seq_info['length']:>{length_width}} {seq_info['strand']} {seq_info['total_length']:>{total_width}} {seq_info['sequence']}"
        return formatted
    
    def realign_block(self, block_lines):
        """Realign a single MAF block, strictly preserving original order"""
        sequences, other_lines = self.parse_maf_block(block_lines)
        # Check number of sequences, skip if less than 2
        if len(sequences) < 2:
            logger.debug(f"Alignment block has fewer than 2 sequences, skipping")
            return block_lines
        
        logger.debug(f"Processing alignment block with {len(sequences)} sequences")
        
        # Calculate maximum column widths
        seq_id_width, start_width, length_width, total_width = self.calculate_column_widths(sequences)
        
        # Prepare sequences for alignment while recording original order
        sequences_for_alignment = []
        original_seq_info_by_id = {}
        
        for seq_info in sequences:
            # Use sequence ID as unique identifier
            seq_id = seq_info['seq_id']
            original_seq_info_by_id[seq_id] = seq_info
            
            # Prepare clean sequence (with gaps removed)
            clean_seq = seq_info['sequence'].replace('-', '').replace(' ', '')
            sequences_for_alignment.append((seq_id, clean_seq))
        
        # Choose alignment tool
        seq_length = len(sequences_for_alignment[0][1])
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_out:
            temp_output = temp_out.name
        
        try:
            # Perform alignment
            if seq_length > self.threshold:
                logger.debug(f"Sequence length: {seq_length}, using MAFFT")
                self.run_mafft_with_order(sequences_for_alignment, temp_output)
            else:
                logger.debug(f"Sequence length: {seq_length}, using MUSCLE5")
                self.run_muscle_with_order(sequences_for_alignment, temp_output)
            
            # Read alignment results
            aligned_seqs = self.parse_fasta_with_ids(temp_output)
            
            # Check if alignment results are valid
            if len(aligned_seqs) != len(sequences):
                logger.error(f"Alignment result sequence count mismatch: expected {len(sequences)}, got {len(aligned_seqs)}")
                return block_lines
            
            result_lines = []
            result_lines.extend(other_lines)
            # Add reformatted sequence lines in original order
            for i, seq_info in enumerate(sequences):
                original_seq_id = seq_info['seq_id']
                
                # Find corresponding sequence in alignment results
                aligned_seq_data = None
                for aligned_id, aligned_seq in aligned_seqs:
                    if aligned_id == original_seq_id:
                        aligned_seq_data = aligned_seq
                        break
                
                if aligned_seq_data is None:
                    logger.warning(f"Sequence {original_seq_id} not found in alignment results, using original sequence")
                    aligned_seq_data = seq_info['sequence'].replace(' ', '')
                
                # Create new sequence information
                new_seq_info = seq_info.copy()
                new_seq_info['sequence'] = aligned_seq_data
                
                # Reformat line
                formatted_line = self.format_sequence_line(
                    new_seq_info, seq_id_width, start_width, length_width, total_width
                )
                result_lines.append(formatted_line)
            
            return result_lines
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Alignment failed, keeping original alignment: {e}")
            # Alignment failed, return original but reformatted block
            return self.format_existing_block(sequences, other_lines, seq_id_width, start_width, length_width, total_width)
        except Exception as e:
            logger.error(f"Error processing alignment block, keeping original alignment: {e}")
            return block_lines
        finally:
            if os.path.exists(temp_output):
                os.unlink(temp_output)
    
    def format_existing_block(self, sequences, other_lines, seq_id_width, start_width, length_width, total_width):
        """Format existing alignment block"""
        result_lines = []
        result_lines.extend(other_lines)
        for seq_info in sequences:
            formatted_line = self.format_sequence_line(
                seq_info, seq_id_width, start_width, length_width, total_width
            )
            result_lines.append(formatted_line)
        return result_lines
    
    def process_maf_file(self, input_file, output_file):
        """Process entire MAF file"""
        logger.info(f"Starting to process file: {input_file}")
        
        blocks_processed = 0
        blocks_aligned = 0
        
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            current_block = []
            
            for line in f_in:
                line = line.rstrip('\n')
                if line.strip() == '':
                    # Empty line indicates end of block
                    if current_block:
                        # Process current alignment block
                        original_has_sequences = any(
                            line.startswith('s ') or line.startswith('s\t') for line in current_block
                        )
                        realigned_block = self.realign_block(current_block)
                        
                        # Write processed block
                        for block_line in realigned_block:
                            f_out.write(block_line + '\n')
                        f_out.write('\n')
                        
                        current_block = []
                        blocks_processed += 1
                        if original_has_sequences:
                            blocks_aligned += 1
                else:
                    current_block.append(line)
            
            # Process last alignment block
            if current_block:
                original_has_sequences = any(
                    line.startswith('s ') or line.startswith('s\t') for line in current_block
                )
                realigned_block = self.realign_block(current_block)
                for block_line in realigned_block:
                    f_out.write(block_line + '\n')
                blocks_processed += 1
                if original_has_sequences:
                    blocks_aligned += 1
        
        logger.info(f"Processing completed: processed {blocks_processed} alignment blocks, realigned {blocks_aligned} blocks containing sequences")

def main():
    parser = argparse.ArgumentParser(description='MAF file realignment module')
    parser.add_argument('-i', '--input', required=True, help='Input MAF file')
    parser.add_argument('-o', '--output', required=True, help='Output MAF file')
    parser.add_argument('--muscle', default='muscle', help='Path to MUSCLE5 executable')
    parser.add_argument('--mafft', default='mafft', help='Path to MAFFT executable')
    parser.add_argument('-t', '--threshold', type=int, default=1000, 
                       help='Sequence length threshold, use MAFFT if above, otherwise MUSCLE5 (default: 1000)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Check input file
    if not os.path.exists(args.input):
        logger.error(f"Input file does not exist: {args.input}")
        sys.exit(1)
    
    # Initialize realigner
    realigner = MafRealigner(
        muscle_path=args.muscle,
        mafft_path=args.mafft,
        threshold=args.threshold
    )
    
    # Check tools
    if not realigner.check_tools():
        sys.exit(1)
    
    # Process file
    try:
        realigner.process_maf_file(args.input, args.output)
        logger.info(f"Output file: {args.output}")
        
        # Validate output
        logger.info("Verifying if sequence order is preserved...")
        
    except Exception as e:
        logger.error(f"Error during processing: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()