#!/usr/bin/env python3
"""
Extract FASTA sequences from genome file based on GFF3 coordinates.
Handles strand orientation correctly (reverse complement for - strand).

Usage:
    python3 extract_sequences_from_gff3.py -g genome.fna.gz -f features.gff3 > sequences.fasta
    python3 extract_sequences_from_gff3.py -g genome.fna.gz < features.gff3 > sequences.fasta
"""

import sys
import argparse
import gzip
from Bio import SeqIO
from Bio.Seq import Seq

def load_genome(genome_file):
    """Load genome sequences into memory for fast access"""
    print(f"Loading genome from {genome_file}...", file=sys.stderr)
    
    try:
        if genome_file.endswith('.gz'):
            handle = gzip.open(genome_file, 'rt')
        else:
            handle = open(genome_file, 'r')
        
        # Load all sequences into dictionary
        genome_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        handle.close()
        
        seq_count = len(genome_dict)
        total_length = sum(len(seq.seq) for seq in genome_dict.values())
        print(f"Loaded {seq_count} sequences ({total_length:,} bp total)", file=sys.stderr)
        
        return genome_dict
        
    except Exception as e:
        print(f"Error loading genome: {e}", file=sys.stderr)
        sys.exit(1)

def parse_gff3_line(line):
    """Parse a single GFF3 line into components"""
    if line.startswith('#') or not line.strip():
        return None
    
    parts = line.strip().split('\t')
    if len(parts) != 9:
        print(f"Warning: Malformed GFF3 line (expected 9 columns): {line}", file=sys.stderr)
        return None
    
    try:
        return {
            'seqid': parts[0],
            'source': parts[1],
            'type': parts[2],
            'start': int(parts[3]),  # GFF3 is 1-based
            'end': int(parts[4]),
            'score': parts[5],
            'strand': parts[6],
            'phase': parts[7],
            'attributes': parts[8]
        }
    except ValueError as e:
        print(f"Warning: Invalid coordinates in GFF3 line: {line}, error: {e}", file=sys.stderr)
        return None

def parse_attributes(attr_string):
    """Parse GFF3 attributes into a dictionary"""
    attributes = {}
    if not attr_string or attr_string == '.':
        return attributes
    
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key] = value
    
    return attributes

def create_fasta_header(feature, seq_length, species_name, feature_counter):
    """Create simple FASTA header: >Species_name-N|type|seqid|start|end|strand"""
    seqid = feature['seqid']
    feature_type = feature['type']
    start = feature['start']
    end = feature['end']
    strand = feature['strand']
    
    # Format: >Drosophila_melanogaster-1|U3|AE014298.5|10348023|10348123|+
    return f"{species_name}-{feature_counter}|{feature_type}|{seqid}|{start}|{end}|{strand}"

def extract_sequence(feature, genome_dict):
    """Extract sequence for a GFF3 feature"""
    seqid = feature['seqid']
    start = feature['start']  # GFF3 is 1-based
    end = feature['end']
    strand = feature['strand']
    
    # Check if sequence exists in genome
    if seqid not in genome_dict:
        print(f"Warning: Sequence '{seqid}' not found in genome, skipping feature", file=sys.stderr)
        return None, None
    
    genome_seq = genome_dict[seqid].seq
    seq_length = len(genome_seq)
    
    # Validate coordinates
    if start < 1 or end > seq_length:
        print(f"Warning: Coordinates {start}-{end} out of bounds for sequence '{seqid}' (length {seq_length}), skipping", file=sys.stderr)
        return None, None
    
    if start > end:
        print(f"Warning: Invalid coordinates {start}-{end} (start > end), skipping", file=sys.stderr)
        return None, None
    
    # Extract sequence (convert to 0-based for Python slicing)
    extracted_seq = genome_seq[start-1:end]
    
    # Handle strand orientation
    if strand == '-':
        extracted_seq = extracted_seq.reverse_complement()
    elif strand != '+' and strand != '.':
        print(f"Warning: Unknown strand '{strand}', treating as forward strand", file=sys.stderr)
    
    return extracted_seq, len(extracted_seq)

def process_gff3(gff3_handle, genome_dict, output_handle, species_name):
    """Process GFF3 file and extract sequences"""
    features_processed = 0
    sequences_extracted = 0
    
    for line_num, line in enumerate(gff3_handle, 1):
        line = line.rstrip('\n\r')
        
        # Skip comments and empty lines
        if line.startswith('#') or not line.strip():
            continue
        
        # Parse feature
        feature = parse_gff3_line(line)
        if feature is None:
            continue
        
        features_processed += 1
        
        # Extract sequence
        sequence, seq_length = extract_sequence(feature, genome_dict)
        
        if sequence is not None:
            # Create FASTA header
            header = create_fasta_header(feature, seq_length, species_name, sequences_extracted + 1)
            
            # Write FASTA entry
            output_handle.write(f">{header}\n")
            
            # Write sequence in lowercase on single line
            seq_str = str(sequence).lower()
            output_handle.write(seq_str + '\n')
            
            sequences_extracted += 1
        
        if features_processed % 1000 == 0:
            print(f"Processed {features_processed} features, extracted {sequences_extracted} sequences...", 
                  file=sys.stderr)
    
    print(f"Total: processed {features_processed} features, extracted {sequences_extracted} sequences", 
          file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Extract FASTA sequences from genome file based on GFF3 coordinates',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract sequences from GFF3 features
  python3 extract_sequences_from_gff3.py -g genome.fna.gz -s Drosophila_melanogaster -f features.gff3 > sequences.fasta
  
  # Using piped GFF3 input
  cat features.gff3 | python3 extract_sequences_from_gff3.py -g genome.fna.gz -s Drosophila_melanogaster > sequences.fasta
  
  # Complete workflow: extract upstream regions and their sequences
  python3 extract_upstream_regions.py -l lengths.tsv < annotations.gff3 | \\
      python3 extract_sequences_from_gff3.py -g genome.fna.gz -s Drosophila_melanogaster > upstream_sequences.fasta
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
    
    # Load genome
    genome_dict = load_genome(args.genome)
    
    # Handle GFF3 input
    if args.gff3:
        gff3_handle = open(args.gff3, 'r')
    else:
        gff3_handle = sys.stdin
    
    # Handle output
    if args.output:
        output_handle = open(args.output, 'w')
    else:
        output_handle = sys.stdout
    
    try:
        process_gff3(gff3_handle, genome_dict, output_handle, args.species)
    finally:
        # Clean up file handles
        if args.gff3:
            gff3_handle.close()
        if args.output:
            output_handle.close()

if __name__ == "__main__":
    main()