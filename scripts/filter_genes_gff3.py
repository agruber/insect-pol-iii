#!/usr/bin/env python3
"""
Filter out specific genes from GFF3 files
Reads GFF3 from stdin and removes specified gene types
"""

import sys
import argparse
from typing import Set

def filter_gff3(genes_to_remove: Set[str]):
    """Filter GFF3 from stdin to stdout, removing specified genes"""
    
    lines_read = 0
    lines_kept = 0
    lines_filtered = 0
    genes_removed_counts = {gene: 0 for gene in genes_to_remove}
    
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
        
        # Check if we have enough fields for a valid GFF3 line
        if len(fields) < 3:
            print(line)  # Keep malformed lines
            lines_kept += 1
            continue
        
        # Check if this gene should be filtered out
        gene_type = fields[2]  # The type field (3rd column) contains the gene name
        if gene_type in genes_to_remove:
            lines_filtered += 1
            genes_removed_counts[gene_type] += 1
        else:
            print(line)
            lines_kept += 1
    
    # Report statistics to stderr
    print(f"Processed {lines_read} lines", file=sys.stderr)
    print(f"Kept {lines_kept} lines", file=sys.stderr)
    print(f"Filtered {lines_filtered} lines", file=sys.stderr)
    
    # Report which genes were removed and how many
    removed_genes = [(gene, count) for gene, count in genes_removed_counts.items() if count > 0]
    if removed_genes:
        print(f"Removed genes:", file=sys.stderr)
        for gene, count in sorted(removed_genes):
            print(f"  {gene}: {count} entries", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Filter out specific genes from GFF3 files (reads from stdin, writes to stdout)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Remove a single gene type
  cat hits.gff3 | python3 filter_genes_gff3.py -r RNaseP_nuc > filtered.gff3
  
  # Remove multiple gene types
  cat hits.gff3 | python3 filter_genes_gff3.py -r RNaseP_nuc -r U6 -r tRNAsec > filtered.gff3
  
  # In a pipeline
  zcat results.tblout.gz | python3 tblout_to_gff3.py | python3 filter_genes_gff3.py -r RNaseP_nuc -r U6 > filtered.gff3
  
  # Combined with e-value filtering
  zcat results.tblout.gz | python3 tblout_to_gff3.py | python3 filter_gff3.py -e 1e-10 | python3 filter_genes_gff3.py -r RNaseP_nuc > filtered.gff3

Note: The gene name is taken from the 3rd column (type) of the GFF3 file
        """
    )
    
    parser.add_argument('-r', '--remove', action='append', dest='genes_to_remove',
                       help='Gene name to filter out (can be used multiple times)')
    
    args = parser.parse_args()
    
    # Check if any genes were specified
    if not args.genes_to_remove:
        print("Error: No genes specified to remove. Use -r option to specify genes.", file=sys.stderr)
        print("Example: python3 filter_genes_gff3.py -r RNaseP_nuc", file=sys.stderr)
        sys.exit(1)
    
    # Convert to set for efficient lookup
    genes_to_remove = set(args.genes_to_remove)
    
    # Run the filter
    try:
        filter_gff3(genes_to_remove)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()