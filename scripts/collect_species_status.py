#!/usr/bin/env python3
"""
Collect status information for all species in the genomes directory.
Outputs a TSV table with status of searches, files, and ncRNAs found.
"""

import sys
import os
import gzip
from pathlib import Path
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import argparse

def check_tblout_file(file_path: Path) -> tuple:
    """Check if tblout file exists and count hits"""
    if not file_path.exists():
        return False, 0, set()

    try:
        hit_count = 0
        ncrna_types = set()
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    hit_count += 1
                    ncrna_types.add(parts[2])  # Target name (ncRNA type)
        return True, hit_count, ncrna_types
    except Exception as e:
        print(f"Error reading {file_path}: {e}", file=sys.stderr)
        return True, 0, set()

def check_blast_file(file_path: Path) -> tuple:
    """Check if blast file exists and count hits"""
    if not file_path.exists():
        return False, 0

    try:
        hit_count = 0
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    hit_count += 1
        return True, hit_count
    except Exception as e:
        print(f"Error reading {file_path}: {e}", file=sys.stderr)
        return True, 0

def check_fasta_file(file_path: Path) -> tuple:
    """Check if FASTA file exists and count sequences"""
    if not file_path.exists():
        return False, 0, set()

    try:
        seq_count = 0
        ncrna_types = set()
        for record in SeqIO.parse(str(file_path), "fasta"):
            seq_count += 1
            # Extract ncRNA type from header (format: species-id|ncRNA_type|...)
            parts = record.id.split('|')
            if len(parts) >= 2:
                ncrna_types.add(parts[1])
        return True, seq_count, ncrna_types
    except Exception as e:
        print(f"Error reading {file_path}: {e}", file=sys.stderr)
        return True, 0, set()

def get_species_status(species_dir: Path) -> dict:
    """Collect all status information for a species"""
    status = {
        'species': species_dir.name,
        'genome_fna': species_dir.joinpath('genome.fna.gz').exists(),
        'genome_tsv': species_dir.joinpath('genome.tsv').exists(),
    }

    # Check tblout files
    snrna_path = species_dir / 'snRNAtRNA.tblout.gz'
    status['snRNAtRNA_exists'], status['snRNAtRNA_hits'], snrna_types = check_tblout_file(snrna_path)

    orcd1_path = species_dir / 'OrCD1.tblout.gz'
    status['OrCD1_exists'], status['OrCD1_hits'], orcd1_types = check_tblout_file(orcd1_path)

    # Check blast file
    blast_path = species_dir / 'noeCR34335.blast.gz'
    status['noeCR34335_exists'], status['noeCR34335_hits'] = check_blast_file(blast_path)

    # Check upstream.fa
    upstream_path = species_dir / 'upstream.fa'
    status['upstream_exists'], status['upstream_seqs'], upstream_types = check_fasta_file(upstream_path)

    # Check annotated.fa
    annotated_path = species_dir / 'annotated.fa'
    status['annotated_exists'], status['annotated_seqs'], _ = check_fasta_file(annotated_path)

    # Check insufficient_diversity.flag
    flag_path = species_dir / 'insufficient_diversity.flag'
    status['insufficient_diversity'] = flag_path.exists()

    # Combine all ncRNA types found
    all_ncrna_types = snrna_types | orcd1_types | upstream_types
    if status['noeCR34335_exists'] and status['noeCR34335_hits'] > 0:
        all_ncrna_types.add('noeCR34335')

    status['ncrna_types'] = ','.join(sorted(all_ncrna_types)) if all_ncrna_types else ''
    status['ncrna_count'] = len(all_ncrna_types)

    return status

def format_status_table(statuses: list) -> str:
    """Format status information as a nice table"""
    df = pd.DataFrame(statuses)

    # Reorder columns for better readability
    column_order = [
        'species',
        'genome_fna', 'genome_tsv',
        'snRNAtRNA_exists', 'snRNAtRNA_hits',
        'OrCD1_exists', 'OrCD1_hits',
        'noeCR34335_exists', 'noeCR34335_hits',
        'upstream_exists', 'upstream_seqs',
        'annotated_exists', 'annotated_seqs',
        'insufficient_diversity',
        'ncrna_count', 'ncrna_types'
    ]

    df = df[column_order]

    # Rename columns for display
    display_names = {
        'species': 'Species',
        'genome_fna': 'Genome',
        'genome_tsv': 'TSV',
        'snRNAtRNA_exists': 'snRNA',
        'snRNAtRNA_hits': '#Hits',
        'OrCD1_exists': 'OrCD1',
        'OrCD1_hits': '#Hits',
        'noeCR34335_exists': 'BLAST',
        'noeCR34335_hits': '#Hits',
        'upstream_exists': 'Upstream',
        'upstream_seqs': '#Seqs',
        'annotated_exists': 'Annotated',
        'annotated_seqs': '#Seqs',
        'insufficient_diversity': 'LowDiv',
        'ncrna_count': '#ncRNA',
        'ncrna_types': 'ncRNA Types'
    }

    df = df.rename(columns=display_names)

    # Convert boolean values to symbols
    bool_cols = ['Genome', 'TSV', 'snRNA', 'OrCD1', 'BLAST', 'Upstream', 'Annotated', 'LowDiv']
    for col in bool_cols:
        if col in df.columns:
            df[col] = df[col].map({True: '✓', False: '-'})

    return df

def main():
    parser = argparse.ArgumentParser(
        description='Collect status information for all species'
    )
    parser.add_argument('-d', '--directory', default='genomes',
                       help='Genomes directory (default: genomes)')
    parser.add_argument('-o', '--output', default=None,
                       help='Output TSV file (default: stdout)')
    parser.add_argument('-f', '--format', choices=['tsv', 'table', 'summary'],
                       default='table',
                       help='Output format (default: table)')
    parser.add_argument('-s', '--species', default=None,
                       help='Filter for specific species (comma-separated)')

    args = parser.parse_args()

    genomes_dir = Path(args.directory)
    if not genomes_dir.exists():
        print(f"Error: Directory {genomes_dir} not found", file=sys.stderr)
        sys.exit(1)

    # Collect status for all species
    statuses = []
    species_dirs = sorted([d for d in genomes_dir.iterdir() if d.is_dir()])

    # Filter species if requested
    if args.species:
        species_filter = set(args.species.split(','))
        species_dirs = [d for d in species_dirs if d.name in species_filter]

    print(f"Collecting status for {len(species_dirs)} species...", file=sys.stderr)

    for species_dir in species_dirs:
        status = get_species_status(species_dir)
        statuses.append(status)

    # Format output
    df = format_status_table(statuses)

    if args.format == 'tsv':
        output = df.to_csv(sep='\t', index=False)
    elif args.format == 'table':
        output = df.to_string(index=False)
    elif args.format == 'summary':
        # Provide summary statistics
        total = len(statuses)
        with_genome = sum(1 for s in statuses if s['genome_fna'])
        with_snrna = sum(1 for s in statuses if s['snRNAtRNA_exists'])
        with_orcd1 = sum(1 for s in statuses if s['OrCD1_exists'])
        with_blast = sum(1 for s in statuses if s['noeCR34335_exists'])
        with_upstream = sum(1 for s in statuses if s['upstream_exists'])
        with_annotated = sum(1 for s in statuses if s['annotated_exists'])
        with_lowdiv = sum(1 for s in statuses if s['insufficient_diversity'])

        output = f"""Status Summary:
Total species: {total}
With genome: {with_genome} ({with_genome/total*100:.1f}%)
With snRNAtRNA search: {with_snrna} ({with_snrna/total*100:.1f}%)
With OrCD1 search: {with_orcd1} ({with_orcd1/total*100:.1f}%)
With BLAST search: {with_blast} ({with_blast/total*100:.1f}%)
With upstream.fa: {with_upstream} ({with_upstream/total*100:.1f}%)
With annotated.fa: {with_annotated} ({with_annotated/total*100:.1f}%)
Insufficient diversity: {with_lowdiv} ({with_lowdiv/total*100:.1f}%)

Average ncRNAs per species: {sum(s['ncrna_count'] for s in statuses)/total:.1f}
Species with >0 ncRNAs: {sum(1 for s in statuses if s['ncrna_count'] > 0)}
"""

    # Write output
    if args.output:
        with open(args.output, 'w') as f:
            f.write(output)
        print(f"Output written to {args.output}", file=sys.stderr)
    else:
        print(output)

if __name__ == "__main__":
    main()