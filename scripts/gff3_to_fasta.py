#!/usr/bin/env python3
"""
Extract FASTA sequences from genome file based on GFF3 coordinates.
Handles strand orientation correctly (reverse complement for - strand).

Usage:
    python3 gff3_to_fasta.py -g genome.fna.gz -f features.gff3 > sequences.fasta
    python3 gff3_to_fasta.py -g genome.fna.gz < features.gff3 > sequences.fasta
"""

import sys
from bioseq_lib import (
    load_fasta_sequences, parse_gff3_file, extract_genomic_sequence,
    write_fasta_sequence, progress_reporter, smart_open,
    create_standard_parser
)

def create_fasta_header(feature, seq_length, species_name, feature_counter):
    """Create simple FASTA header: >Species_name-N|type|seqid|start|end|strand"""
    return f"{species_name}-{feature_counter}|{feature.type}|{feature.seqid}|{feature.start}|{feature.end}|{feature.strand}"

def process_gff3(gff3_filename, genome_dict, output_handle, species_name):
    """Process GFF3 file and extract sequences using library functions"""
    features_processed = 0
    sequences_extracted = 0
    
    for feature in parse_gff3_file(gff3_filename):
        features_processed += 1
        
        # Extract sequence using library function
        sequence = extract_genomic_sequence(feature, genome_dict)
        
        if sequence is not None:
            # Create FASTA header
            header = create_fasta_header(feature, len(sequence), species_name, sequences_extracted + 1)
            
            # Write sequence using library function
            write_fasta_sequence(output_handle, header, sequence, lowercase=True)
            sequences_extracted += 1
        
        progress_reporter(features_processed)
    
    print(f"Total: processed {features_processed} features, extracted {sequences_extracted} sequences", 
          file=sys.stderr)

def main():
    parser = create_standard_parser(
        description='Extract FASTA sequences from genome file based on GFF3 coordinates',
        epilog="""
Examples:
  # Extract sequences from GFF3 features
  python3 gff3_to_fasta.py -g genome.fna.gz -s Drosophila_melanogaster -f features.gff3 > sequences.fasta
  
  # Using piped GFF3 input
  cat features.gff3 | python3 gff3_to_fasta.py -g genome.fna.gz -s Drosophila_melanogaster > sequences.fasta
        """
    )
    
    parser.add_argument('-g', '--genome', required=True,
                       help='Input genome FASTA file (supports .gz compression)')
    parser.add_argument('-s', '--species', required=True,
                       help='Species name for FASTA headers (e.g., Drosophila_melanogaster)')
    parser.add_argument('-f', '--gff3',
                       help='Input GFF3 file (use stdin if not provided)')
    parser.add_argument('-o', '--output',
                       help='Output FASTA file (use stdout if not provided)')
    
    args = parser.parse_args()
    
    # Load genome using library function
    genome_dict = load_fasta_sequences(args.genome)
    
    # Handle output
    output_handle = smart_open(args.output, 'w')
    
    try:
        process_gff3(args.gff3, genome_dict, output_handle, args.species)
    finally:
        if args.output:
            output_handle.close()

if __name__ == "__main__":
    main()