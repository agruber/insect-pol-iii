#!/usr/bin/env python3
"""
Extract upstream regions from GFF3 annotations.
Generates new GFF3 entries for regions upstream of existing features.

Usage:
    python3 extract_upstream_regions.py -l lengths.tsv < input.gff3 > upstream.gff3
    python3 extract_upstream_regions.py -l lengths.tsv -u 200 input.gff3 > upstream.gff3
"""

import sys
from bioseq_lib import (
    smart_open, parse_gff3_file, calculate_upstream_region,
    progress_reporter, create_standard_parser
)

def load_sequence_lengths(lengths_file):
    """Load sequence lengths from TSV file"""
    lengths = {}
    try:
        with smart_open(lengths_file, 'r') as f:
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
        
    except Exception as e:
        print(f"Error reading lengths file: {e}", file=sys.stderr)
        sys.exit(1)

def process_gff3(input_filename, output_handle, seq_lengths, upstream_distance):
    """Process GFF3 input and generate upstream regions using library functions"""
    features_processed = 0
    upstream_created = 0
    
    # Write GFF3 header
    output_handle.write("##gff-version 3\n")
    output_handle.write(f"##upstream-regions-distance {upstream_distance}\n")
    
    for feature in parse_gff3_file(input_filename):
        features_processed += 1
        
        # Create upstream region using library function
        upstream_feature = calculate_upstream_region(feature, upstream_distance, seq_lengths)
        
        if upstream_feature:
            # Write upstream region
            output_handle.write(upstream_feature.to_line() + '\n')
            upstream_created += 1
        
        progress_reporter(features_processed)
    
    print(f"Total: processed {features_processed} features, created {upstream_created} upstream regions", 
          file=sys.stderr)

def main():
    parser = create_standard_parser(
        description='Extract upstream regions from GFF3 annotations',
        epilog="""
Examples:
  # Basic usage with piped input
  python3 extract_upstream_regions.py -l lengths.tsv < input.gff3 > upstream.gff3
  
  # With custom upstream distance
  python3 extract_upstream_regions.py -l lengths.tsv -u 200 input.gff3 > upstream.gff3
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
    
    # Handle output
    output_handle = smart_open(args.output, 'w')
    
    try:
        process_gff3(args.input_gff3, output_handle, seq_lengths, args.upstream)
    finally:
        if args.output:
            output_handle.close()

if __name__ == "__main__":
    main()