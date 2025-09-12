#!/usr/bin/env python3
"""
Filter GFF3 files by e-value thresholds
Reads GFF3 from stdin and filters based on e-value in the score column
Supports global e-value cutoff and gene-specific cutoffs
"""

import sys
import argparse
from typing import Dict, Optional

def parse_gene_filters(filter_strings):
    """Parse gene-specific e-value filters from command line arguments"""
    gene_filters = {}
    for filter_str in filter_strings:
        parts = filter_str.split(':')
        if len(parts) != 2:
            raise ValueError(f"Invalid gene filter format: {filter_str}. Expected format: gene_name:evalue")
        gene_name = parts[0]
        try:
            evalue = float(parts[1])
        except ValueError:
            raise ValueError(f"Invalid e-value in filter: {filter_str}")
        gene_filters[gene_name] = evalue
    return gene_filters

def should_keep_line(fields, default_evalue: float, gene_filters: Dict[str, float]) -> bool:
    """Determine if a GFF3 line should be kept based on e-value filters"""
    if len(fields) < 6:
        return True  # Keep malformed lines or headers
    
    gene_type = fields[2]  # The type field (3rd column) contains the gene name
    try:
        evalue = float(fields[5])  # The score field (6th column) contains the e-value
    except (ValueError, IndexError):
        return True  # Keep lines where e-value cannot be parsed
    
    # Check gene-specific filter first
    if gene_type in gene_filters:
        threshold = gene_filters[gene_type]
    else:
        threshold = default_evalue
    
    return evalue <= threshold

def filter_gff3(default_evalue: float = 1e-5, gene_filters: Optional[Dict[str, float]] = None):
    """Filter GFF3 from stdin to stdout based on e-value thresholds"""
    
    if gene_filters is None:
        gene_filters = {}
    
    lines_read = 0
    lines_kept = 0
    lines_filtered = 0
    
    for line in sys.stdin:
        line = line.strip()
        lines_read += 1
        
        # Always keep empty lines and comments/headers
        if not line or line.startswith('#'):
            print(line)
            lines_kept += 1
            continue
        
        # Parse GFF3 line
        fields = line.split('\t')
        
        # Apply filter
        if should_keep_line(fields, default_evalue, gene_filters):
            print(line)
            lines_kept += 1
        else:
            lines_filtered += 1
    
    # Report statistics to stderr
    print(f"Processed {lines_read} lines", file=sys.stderr)
    print(f"Kept {lines_kept} lines", file=sys.stderr)
    print(f"Filtered {lines_filtered} lines", file=sys.stderr)
    
    if gene_filters:
        print(f"Applied gene-specific filters for: {', '.join(gene_filters.keys())}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Filter GFF3 files by e-value thresholds (reads from stdin, writes to stdout)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Filter with global e-value cutoff
  cat hits.gff3 | python3 filter_gff3.py -e 1e-10 > filtered.gff3
  
  # Filter with gene-specific e-value cutoffs
  cat hits.gff3 | python3 filter_gff3.py -e 1e-5 -g U6:1e-20 -g U1:1e-15 > filtered.gff3
  
  # Multiple gene-specific filters
  cat hits.gff3 | python3 filter_gff3.py -e 1e-3 -g U6:1e-20 -g U1:1e-15 -g U2:1e-18 > filtered.gff3
  
  # In a pipeline
  zcat results.tblout.gz | python3 tblout_to_gff3.py | python3 filter_gff3.py -e 1e-10 > filtered.gff3

Note: The gene name is taken from the 3rd column (type) of the GFF3 file
      The e-value is expected to be in the 6th column (score) of the GFF3 file
        """
    )
    
    parser.add_argument('-e', '--evalue', type=float, default=1e-5,
                       help='Default e-value cutoff for filtering (default: 1e-5)')
    parser.add_argument('-g', '--gene', action='append', dest='gene_filters',
                       help='Gene-specific e-value filter in format gene_name:evalue (can be used multiple times)')
    
    args = parser.parse_args()
    
    # Parse gene-specific filters
    gene_filters = {}
    if args.gene_filters:
        try:
            gene_filters = parse_gene_filters(args.gene_filters)
        except ValueError as e:
            print(f"Error parsing gene filters: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Run the filter
    try:
        filter_gff3(args.evalue, gene_filters)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()