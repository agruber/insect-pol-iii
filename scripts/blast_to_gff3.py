#!/usr/bin/env python3
"""
Convert BLAST tabular output to GFF3 format
Reads BLAST format 6 output from stdin and converts to GFF3
"""

import sys
import argparse
from datetime import datetime

def parse_blast_line(line):
    """Parse a BLAST tabular output line into components"""
    fields = line.strip().split('\t')
    if len(fields) < 15:
        return None
    
    try:
        return {
            'qseqid': fields[0],      # Query sequence ID
            'sseqid': fields[1],      # Subject sequence ID  
            'pident': float(fields[2]), # Percent identity
            'length': int(fields[3]),   # Alignment length
            'mismatch': int(fields[4]), # Number of mismatches
            'gapopen': int(fields[5]),  # Number of gap openings
            'qstart': int(fields[6]),   # Query start
            'qend': int(fields[7]),     # Query end
            'sstart': int(fields[8]),   # Subject start
            'send': int(fields[9]),     # Subject end
            'evalue': float(fields[10]), # E-value
            'bitscore': float(fields[11]), # Bit score
            'qlen': int(fields[12]),    # Query length
            'slen': int(fields[13]),    # Subject length
            'sstrand': fields[14] if len(fields) > 14 else 'plus'  # Subject strand
        }
    except (ValueError, IndexError):
        return None

def blast_to_gff3_entry(hit, source="blast", gene_name=None):
    """Convert a BLAST hit to GFF3 format"""
    
    # Handle coordinates and strand
    if hit['sstart'] <= hit['send']:
        start = hit['sstart']
        end = hit['send']
        strand = '+' if hit['sstrand'] == 'plus' else '+'
    else:
        start = hit['send']
        end = hit['sstart']
        strand = '-' if hit['sstrand'] == 'minus' else '-'
    
    # Handle strand information from sstrand field
    if hit['sstrand'] == 'minus':
        strand = '-'
    elif hit['sstrand'] == 'plus':
        strand = '+'
    else:
        strand = '.'
    
    # Use provided gene name or fall back to query sequence ID
    feature_type = gene_name if gene_name else hit['qseqid']
    
    # Create GFF3 line
    gff3_line = '\t'.join([
        hit['sseqid'],               # seqid (subject sequence)
        source,                      # source
        feature_type,               # type (gene name or query sequence name)
        str(start),                 # start (1-based)
        str(end),                   # end (1-based)
        str(hit['evalue']),         # score (using e-value)
        strand,                     # strand
        '.',                        # phase (not applicable)
        '.'                         # attributes (empty)
    ])
    
    return gff3_line

def convert_blast_to_gff3(source="blast", gene_name=None):
    """Convert BLAST results from stdin to GFF3 format, writing to stdout"""
    
    # Read from stdin
    in_fh = sys.stdin
    
    # Always write to stdout
    out_fh = sys.stdout
    
    # Write GFF3 header
    out_fh.write("##gff-version 3\n")
    out_fh.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
    out_fh.write(f"##source-version blast_to_gff3.py\n")
    
    hits_written = 0
    
    for line in in_fh:
        line = line.strip()
        
        # Skip comments and empty lines
        if line.startswith('#') or not line:
            continue
        
        # Parse the hit
        hit = parse_blast_line(line)
        if hit is None:
            continue
        
        # Convert to GFF3 and write
        gff3_line = blast_to_gff3_entry(hit, source, gene_name)
        out_fh.write(gff3_line + '\n')
        hits_written += 1
    
    print(f"Converted {hits_written} BLAST hits to GFF3 format", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Convert BLAST tabular output to GFF3 format (reads from stdin, outputs to stdout)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert BLAST output to GFF3
  cat results.blast | python3 blast_to_gff3.py > hits.gff3
  
  # From compressed file
  zcat results.blast.gz | python3 blast_to_gff3.py > hits.gff3
  
  # In a pipeline with processing
  zcat results.blast.gz | python3 process_blast_hits.py | python3 blast_to_gff3.py > filtered_hits.gff3
  
  # With custom source field
  zcat results.blast.gz | python3 blast_to_gff3.py -s blastn > hits.gff3

Input format expected (BLAST format 6 with extended fields):
  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand

Output:
  Standard GFF3 format with:
  - seqid: subject sequence ID
  - source: blast (or custom)
  - type: query sequence ID
  - coordinates: subject start/end
  - score: e-value
  - strand: derived from subject strand
  - attributes: empty (.)
        """
    )
    
    parser.add_argument('-s', '--source', default='blast',
                       help='Source field for GFF3 (default: blast)')
    parser.add_argument('-g', '--gene-name',
                       help='Gene name to use in type field (overrides query sequence ID)')
    
    args = parser.parse_args()
    
    try:
        convert_blast_to_gff3(args.source, args.gene_name)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()