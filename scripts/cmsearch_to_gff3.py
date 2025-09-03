#!/usr/bin/env python3
"""
Convert cmsearch results to GFF3 format
Handles coordinate system conversion and E-value filtering
"""

import sys
import argparse
import gzip
from datetime import datetime

def parse_cmsearch_line(line):
    """Parse a cmsearch result line into components"""
    fields = line.strip().split()
    if len(fields) < 16:
        return None
    
    return {
        'target_name': fields[0],
        'target_accession': fields[1] if fields[1] != '-' else None,
        'query_name': fields[2],
        'query_accession': fields[3],
        'model_type': fields[4],
        'model_from': int(fields[5]),
        'model_to': int(fields[6]),
        'seq_from': int(fields[7]),
        'seq_to': int(fields[8]),
        'strand': fields[9],
        'trunc': fields[10],
        'pass': int(fields[11]),
        'gc': float(fields[12]),
        'bias': float(fields[13]),
        'score': float(fields[14]),
        'evalue': float(fields[15]),
        'inc': fields[16],
        'description': ' '.join(fields[17:]) if len(fields) > 17 else ''
    }

def cmsearch_to_gff3_entry(hit, source="cmsearch"):
    """Convert a cmsearch hit to GFF3 format"""
    
    # Handle coordinates: cmsearch uses 1-based coordinates
    # GFF3 uses 1-based coordinates, so no conversion needed for positions
    # But we need to ensure start <= end for GFF3
    if hit['strand'] == '+':
        start = hit['seq_from']
        end = hit['seq_to']
    else:  # strand == '-'
        start = hit['seq_to']
        end = hit['seq_from']
    
    # Ensure start <= end (GFF3 requirement)
    if start > end:
        start, end = end, start
    
    # Build attributes
    attributes = [
        f"ID={hit['query_name']}_{hit['target_name']}_{start}_{end}",
        f"Name={hit['query_name']}",
        f"Target={hit['query_name']} {hit['model_from']} {hit['model_to']}",
        f"score={hit['score']}",
        f"evalue={hit['evalue']}",
        f"bias={hit['bias']}",
        f"gc_content={hit['gc']}"
    ]
    
    if hit['query_accession'] != '-':
        attributes.append(f"query_accession={hit['query_accession']}")
    
    if hit['description']:
        attributes.append(f"Note={hit['description']}")
    
    attributes_str = ';'.join(attributes)
    
    # Create GFF3 line
    gff3_line = '\t'.join([
        hit['target_name'],          # seqid
        source,                      # source  
        hit['query_name'],           # type
        str(start),                  # start (1-based)
        str(end),                    # end (1-based)
        str(hit['score']),           # score
        hit['strand'],               # strand
        '.',                         # phase (not applicable for ncRNA)
        attributes_str               # attributes
    ])
    
    return gff3_line

def convert_cmsearch_to_gff3(input_file, output_file, evalue_cutoff=1e-5, source="cmsearch"):
    """Convert cmsearch results file to GFF3 format"""
    
    # Determine if input is gzipped
    if input_file.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    # Determine output file handle
    if output_file == '-':
        out_fh = sys.stdout
    else:
        out_fh = open(output_file, 'w')
    
    try:
        with open_func(input_file, mode) as in_fh:
            # Write GFF3 header
            out_fh.write("##gff-version 3\n")
            out_fh.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
            out_fh.write(f"##source-version cmsearch_to_gff3.py\n")
            out_fh.write(f"##evalue-cutoff {evalue_cutoff}\n")
            
            hits_written = 0
            hits_filtered = 0
            
            for line in in_fh:
                line = line.strip()
                
                # Skip comments and empty lines
                if line.startswith('#') or not line:
                    continue
                
                # Parse the hit
                hit = parse_cmsearch_line(line)
                if hit is None:
                    continue
                
                # Apply E-value filter
                if hit['evalue'] > evalue_cutoff:
                    hits_filtered += 1
                    continue
                
                # Convert to GFF3 and write
                gff3_line = cmsearch_to_gff3_entry(hit, source)
                out_fh.write(gff3_line + '\n')
                hits_written += 1
            
    finally:
        if output_file != '-':
            out_fh.close()
    
    print(f"Converted {hits_written} hits to GFF3 format", file=sys.stderr)
    print(f"Filtered {hits_filtered} hits with E-value > {evalue_cutoff}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Convert cmsearch results to GFF3 format')
    parser.add_argument('input_file', help='Input cmsearch results file (can be .gz)')
    parser.add_argument('-o', '--output', default='-', 
                       help='Output GFF3 file (default: stdout)')
    parser.add_argument('-e', '--evalue', type=float, default=1e-5,
                       help='E-value cutoff for filtering hits (default: 1e-5)')
    parser.add_argument('-s', '--source', default='cmsearch',
                       help='Source field for GFF3 (default: cmsearch)')
    
    args = parser.parse_args()
    
    try:
        convert_cmsearch_to_gff3(args.input_file, args.output, args.evalue, args.source)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()