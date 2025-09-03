#!/usr/bin/env python3
"""
Extract sequence lengths from FASTA files.
Supports both file input and stdin (for piped gzipped files).

Usage:
    zcat genome.fna.gz | python3 get_sequence_lengths.py > lengths.tsv
    python3 get_sequence_lengths.py genome.fna > lengths.tsv
    python3 get_sequence_lengths.py -o lengths.tsv genome.fna

Output format: sequence_id<tab>length
"""

import sys
from bioseq_lib import (
    smart_open, progress_reporter, parse_fasta_from_handle,
    create_standard_parser, add_io_arguments
)

def get_sequence_lengths(input_handle, output_handle):
    """Extract sequence lengths from FASTA input"""
    try:
        sequences = parse_fasta_from_handle(input_handle)
        count = 0
        
        for record in sequences:
            seq_id = record.id
            length = len(record.seq)
            output_handle.write(f"{seq_id}\t{length}\n")
            count += 1
            
            progress_reporter(count)
        
        print(f"Total sequences processed: {count}", file=sys.stderr)
        
    except Exception as e:
        print(f"Error processing FASTA: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = create_standard_parser(
        description='Extract sequence lengths from FASTA files',
        epilog="""
Examples:
  # From gzipped file via pipe
  zcat genome.fna.gz | python3 get_sequence_lengths.py > lengths.tsv
  
  # From regular file
  python3 get_sequence_lengths.py genome.fna > lengths.tsv
  
  # With output file specified
  python3 get_sequence_lengths.py -o lengths.tsv genome.fna
        """
    )
    
    add_io_arguments(parser, 
                    'Input FASTA file', 
                    'Output TSV file')
    
    args = parser.parse_args()
    
    # Handle I/O using library functions
    input_handle = smart_open(args.input, 'r')
    output_handle = smart_open(args.output, 'w')
    
    try:
        get_sequence_lengths(input_handle, output_handle)
    finally:
        if args.input:
            input_handle.close()
        if args.output:
            output_handle.close()

if __name__ == "__main__":
    main()