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
import argparse
import gzip
from Bio import SeqIO

def get_sequence_lengths(input_handle, output_handle):
    """Extract sequence lengths from FASTA input using BioPython"""
    try:
        sequences = SeqIO.parse(input_handle, 'fasta')
        count = 0
        
        for record in sequences:
            seq_id = record.id
            length = len(record.seq)
            output_handle.write(f"{seq_id}\t{length}\n")
            count += 1
            
            if count % 1000 == 0:
                print(f"Processed {count} sequences...", file=sys.stderr)
        
        print(f"Total sequences processed: {count}", file=sys.stderr)
        
    except Exception as e:
        print(f"Error processing FASTA: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Extract sequence lengths from FASTA files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
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
    
    parser.add_argument('input_file', nargs='?', 
                       help='Input FASTA file (use stdin if not provided)')
    parser.add_argument('-o', '--output', 
                       help='Output file (use stdout if not provided)')
    
    args = parser.parse_args()
    
    # Handle input
    if args.input_file:
        if args.input_file.endswith('.gz'):
            input_handle = gzip.open(args.input_file, 'rt')
        else:
            input_handle = open(args.input_file, 'r')
    else:
        input_handle = sys.stdin
    
    # Handle output
    if args.output:
        output_handle = open(args.output, 'w')
    else:
        output_handle = sys.stdout
    
    try:
        get_sequence_lengths(input_handle, output_handle)
    finally:
        # Clean up file handles
        if args.input_file:
            input_handle.close()
        if args.output:
            output_handle.close()

if __name__ == "__main__":
    main()