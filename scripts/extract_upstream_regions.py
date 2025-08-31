#!/usr/bin/env python3
"""
Extract upstream regions from GFF3 annotations.
Generates new GFF3 entries for regions upstream of existing features.

Usage:
    python3 extract_upstream_regions.py -l lengths.tsv < input.gff3 > upstream.gff3
    python3 extract_upstream_regions.py -l lengths.tsv -u 200 input.gff3 > upstream.gff3

The script handles strand orientation correctly:
- For + strand: upstream = [start - upstream_distance, start - 1]  
- For - strand: upstream = [end + 1, end + upstream_distance]
"""

import sys
import argparse
import re

def load_sequence_lengths(lengths_file):
    """Load sequence lengths from TSV file"""
    lengths = {}
    try:
        with open(lengths_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) != 2:
                    print(f"Warning: Malformed line {line_num} in {lengths_file}: {line}", 
                          file=sys.stderr)
                    continue
                
                seq_id, length_str = parts
                try:
                    length = int(length_str)
                    lengths[seq_id] = length
                except ValueError:
                    print(f"Warning: Invalid length '{length_str}' for sequence '{seq_id}' at line {line_num}", 
                          file=sys.stderr)
        
        print(f"Loaded lengths for {len(lengths)} sequences", file=sys.stderr)
        return lengths
        
    except FileNotFoundError:
        print(f"Error: Lengths file '{lengths_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading lengths file: {e}", file=sys.stderr)
        sys.exit(1)

def parse_gff3_line(line):
    """Parse a single GFF3 line into components"""
    if line.startswith('#') or not line.strip():
        return None
    
    parts = line.strip().split('\t')
    if len(parts) != 9:
        print(f"Warning: Malformed GFF3 line (expected 9 columns): {line}", file=sys.stderr)
        return None
    
    return {
        'seqid': parts[0],
        'source': parts[1],
        'type': parts[2],
        'start': int(parts[3]),
        'end': int(parts[4]),
        'score': parts[5],
        'strand': parts[6],
        'phase': parts[7],
        'attributes': parts[8]
    }

def create_upstream_region(feature, upstream_distance, seq_lengths):
    """Create upstream region for a feature"""
    seqid = feature['seqid']
    strand = feature['strand']
    
    # Get sequence length
    if seqid not in seq_lengths:
        print(f"Warning: No length information for sequence '{seqid}', skipping", file=sys.stderr)
        return None
    
    seq_length = seq_lengths[seqid]
    
    # Calculate upstream region coordinates
    if strand == '+':
        # For + strand: upstream is before the start
        upstream_end = feature['start'] - 1
        upstream_start = max(1, upstream_end - upstream_distance + 1)
    elif strand == '-':
        # For - strand: upstream is after the end
        upstream_start = feature['end'] + 1
        upstream_end = min(seq_length, upstream_start + upstream_distance - 1)
    else:
        # Unknown strand - skip
        print(f"Warning: Unknown strand '{strand}' for feature, skipping upstream region", 
              file=sys.stderr)
        return None
    
    # Check if region is valid
    if upstream_start > upstream_end or upstream_start < 1 or upstream_end > seq_length:
        print(f"Warning: Invalid upstream region [{upstream_start}, {upstream_end}] for sequence {seqid} (length {seq_length}), skipping", 
              file=sys.stderr)
        return None
    
    # Create new GFF3 entry for upstream region
    # Modify attributes to indicate this is an upstream region
    new_attributes = feature['attributes']
    
    # Add upstream information to attributes
    if new_attributes and not new_attributes.endswith(';'):
        new_attributes += ';'
    new_attributes += f"upstream_of={feature['type']};upstream_distance={upstream_distance}"
    
    # Create upstream feature with original type name
    upstream_feature = {
        'seqid': seqid,
        'source': feature['source'],
        'type': feature['type'],  # Keep original type name
        'start': upstream_start,
        'end': upstream_end,
        'score': '.',
        'strand': strand,
        'phase': '.',
        'attributes': new_attributes
    }
    
    return upstream_feature

def format_gff3_line(feature):
    """Format a feature dict back to GFF3 line"""
    return f"{feature['seqid']}\t{feature['source']}\t{feature['type']}\t{feature['start']}\t{feature['end']}\t{feature['score']}\t{feature['strand']}\t{feature['phase']}\t{feature['attributes']}"

def process_gff3(input_handle, output_handle, seq_lengths, upstream_distance):
    """Process GFF3 input and generate upstream regions"""
    features_processed = 0
    upstream_created = 0
    
    # Write GFF3 header
    output_handle.write("##gff-version 3\n")
    output_handle.write(f"##upstream-regions-distance {upstream_distance}\n")
    
    for line_num, line in enumerate(input_handle, 1):
        line = line.rstrip('\n\r')
        
        # Pass through comments and empty lines
        if line.startswith('#') or not line.strip():
            output_handle.write(line + '\n')
            continue
        
        # Parse feature line
        feature = parse_gff3_line(line)
        if feature is None:
            continue
        
        features_processed += 1
        
        # Create upstream region
        upstream_feature = create_upstream_region(feature, upstream_distance, seq_lengths)
        
        if upstream_feature:
            # Write only the upstream region
            output_handle.write(format_gff3_line(upstream_feature) + '\n')
            upstream_created += 1
        
        if features_processed % 1000 == 0:
            print(f"Processed {features_processed} features, created {upstream_created} upstream regions...", 
                  file=sys.stderr)
    
    print(f"Total: processed {features_processed} features, created {upstream_created} upstream regions", 
          file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Extract upstream regions from GFF3 annotations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with piped input
  python3 extract_upstream_regions.py -l lengths.tsv < input.gff3 > upstream.gff3
  
  # With custom upstream distance
  python3 extract_upstream_regions.py -l lengths.tsv -u 200 input.gff3 > upstream.gff3
  
  # Generate lengths file first, then extract upstream
  zcat genome.fna.gz | python3 get_sequence_lengths.py > lengths.tsv
  python3 extract_upstream_regions.py -l lengths.tsv < annotations.gff3 > upstream.gff3
        """
    )
    
    parser.add_argument('-l', '--lengths', required=True,
                       help='TSV file with sequence lengths (seqid<tab>length)')
    parser.add_argument('-u', '--upstream', type=int, default=100,
                       help='Upstream distance in nucleotides (default: 100)')
    parser.add_argument('input_gff3', nargs='?',
                       help='Input GFF3 file (use stdin if not provided)')
    parser.add_argument('-o', '--output',
                       help='Output GFF3 file (use stdout if not provided)')
    
    args = parser.parse_args()
    
    if args.upstream <= 0:
        print("Error: Upstream distance must be positive", file=sys.stderr)
        sys.exit(1)
    
    # Load sequence lengths
    seq_lengths = load_sequence_lengths(args.lengths)
    
    # Handle input
    if args.input_gff3:
        input_handle = open(args.input_gff3, 'r')
    else:
        input_handle = sys.stdin
    
    # Handle output  
    if args.output:
        output_handle = open(args.output, 'w')
    else:
        output_handle = sys.stdout
    
    try:
        process_gff3(input_handle, output_handle, seq_lengths, args.upstream)
    finally:
        # Clean up file handles
        if args.input_gff3:
            input_handle.close()
        if args.output:
            output_handle.close()

if __name__ == "__main__":
    main()