#!/usr/bin/env python3
"""
Export all results from results/ folder into a comprehensive TSV file.

Usage:
    python3 scripts/export_results_to_tsv.py -o results_summary.tsv
    python3 scripts/export_results_to_tsv.py --results-dir results --taxonomy data/genomes_annotated.tsv -o output.tsv
"""

import sys
import csv
from pathlib import Path
import argparse
from collections import defaultdict


def load_taxonomy(tsv_file="data/genomes_annotated.tsv"):
    """Load taxonomy information from genomes_annotated.tsv"""
    taxonomy = {}

    tsv_path = Path(tsv_file)
    if not tsv_path.exists():
        print(f"Warning: Taxonomy file not found: {tsv_file}", file=sys.stderr)
        return taxonomy

    with open(tsv_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue

            species = parts[0]
            assembly_url = parts[2] if len(parts) > 2 else ''

            # Extract assembly name from URL (last part after /)
            assembly_name = ''
            if assembly_url:
                assembly_name = assembly_url.rstrip('/').split('/')[-1]

            order = parts[5] if len(parts) > 5 and parts[5] != 'NA' else ''
            suborder = parts[6] if len(parts) > 6 and parts[6] != 'NA' else ''
            infraorder = parts[7] if len(parts) > 7 and parts[7] != 'NA' else ''
            superfamily = parts[8] if len(parts) > 8 and parts[8] != 'NA' else ''
            family = parts[9] if len(parts) > 9 and parts[9] != 'NA' else ''
            genus = parts[10] if len(parts) > 10 and parts[10] != 'NA' else ''
            subgenus = parts[11] if len(parts) > 11 and parts[11] != 'NA' else ''
            group = parts[12] if len(parts) > 12 and parts[12] != 'NA' else ''
            subgroup = parts[13] if len(parts) > 13 and parts[13] != 'NA' else ''

            taxonomy[species] = {
                'order': order,
                'suborder': suborder,
                'infraorder': infraorder,
                'superfamily': superfamily,
                'family': family,
                'genus': genus,
                'subgenus': subgenus,
                'group': group,
                'subgroup': subgroup,
                'assembly_name': assembly_name,
                # Store sort key like concat_pdfs_by_taxonomy.py
                'sort_key': [order, suborder, infraorder, superfamily, family, genus, subgenus, group, subgroup, species]
            }

    return taxonomy


def parse_fasta(fasta_file):
    """Parse FASTA file and return dict of sequences"""
    sequences = {}

    if not fasta_file.exists():
        return sequences

    with open(fasta_file, 'r') as f:
        current_header = None
        current_seq = []

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_header:
                    sequences[current_header] = ''.join(current_seq)

                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget last sequence
        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences


def parse_stats_file(stats_file):
    """Parse lncrna.stats file"""
    stats = {}

    if not stats_file.exists():
        return stats

    with open(stats_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq_name = row['sequence_name']
            stats[seq_name] = {
                'motif_5prime_type': row.get('motif_5prime_type', ''),
                'motif_3prime_type': row.get('motif_3prime_type', ''),
                'longest_polyT_stretch': row.get('longest_polyT_stretch', ''),
                'trailing_T_count': row.get('trailing_T_count', '')
            }

    return stats


def parse_scores_file(scores_file):
    """Parse noe.scores file"""
    scores = {}

    if not scores_file.exists():
        return scores

    with open(scores_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq_id = row['Sequence_ID']
            scores[seq_id] = row.get('E_value', '')

    return scores


def parse_upstream_scores(upstream_file):
    """Parse upstream.score file for Pol III scores"""
    pol3_scores = {}

    if not upstream_file.exists():
        return pol3_scores

    with open(upstream_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq_id = row['ID']
            pol3_scores[seq_id] = row.get('Pol3_Score', '')

    return pol3_scores


def export_results(results_dir, taxonomy, output_file):
    """Export all results to TSV"""

    results_path = Path(results_dir)
    if not results_path.exists():
        print(f"Error: Results directory not found: {results_dir}", file=sys.stderr)
        sys.exit(1)

    # First, collect all data rows with their taxonomy for sorting
    all_rows = []
    total_sequences = 0
    species_count = 0

    # Process each species directory
    for species_dir in results_path.iterdir():
        if not species_dir.is_dir():
            continue

        species_name = species_dir.name.replace('_', ' ')

        # Get taxonomy info
        tax_info = taxonomy.get(species_name, {})

        # Parse all data files
        lncrna_file = species_dir / 'lncrna.fa'
        stats_file = species_dir / 'lncrna.stats'
        scores_file = species_dir / 'noe.scores'
        upstream_file = species_dir / 'upstream.score'

        if not lncrna_file.exists():
            continue

        species_count += 1
        sequences = parse_fasta(lncrna_file)
        stats = parse_stats_file(stats_file)
        noe_scores = parse_scores_file(scores_file)
        pol3_scores = parse_upstream_scores(upstream_file)

        # Process each sequence
        for header, sequence in sequences.items():
            # Parse header: ID|type|contig|start|end|strand
            parts = header.split('|')
            if len(parts) < 6:
                print(f"Warning: Malformed header: {header}", file=sys.stderr)
                continue

            seq_id = parts[0]
            contig = parts[2]
            start = parts[3]
            end = parts[4]
            strand = parts[5]
            length = len(sequence)

            # Get stats for this sequence
            seq_stats = stats.get(header, {})

            # Build row
            row = {
                'order': tax_info.get('order', ''),
                'suborder': tax_info.get('suborder', ''),
                'infraorder': tax_info.get('infraorder', ''),
                'superfamily': tax_info.get('superfamily', ''),
                'family': tax_info.get('family', ''),
                'genus': tax_info.get('genus', ''),
                'subgenus': tax_info.get('subgenus', ''),
                'group': tax_info.get('group', ''),
                'subgroup': tax_info.get('subgroup', ''),
                'Species': species_name,
                'assembly_name': tax_info.get('assembly_name', ''),
                'Sequence_ID': seq_id,
                'Contig': contig,
                'Start': start,
                'End': end,
                'Strand': strand,
                'Length': length,
                'motif_5prime_type': seq_stats.get('motif_5prime_type', ''),
                'motif_3prime_type': seq_stats.get('motif_3prime_type', ''),
                'longest_polyT_stretch': seq_stats.get('longest_polyT_stretch', ''),
                'trailing_T_count': seq_stats.get('trailing_T_count', ''),
                'noe_evalue': noe_scores.get(header, ''),
                'pol3_score': pol3_scores.get(seq_id, ''),
                'Sequence': sequence
            }

            # Store row with sort key
            sort_key = tax_info.get('sort_key', [''] * 9 + [species_name])
            all_rows.append((sort_key, row))
            total_sequences += 1

    # Sort all rows by taxonomy
    all_rows.sort(key=lambda x: x[0])

    # Write sorted data to file
    fieldnames = [
        'order', 'suborder', 'infraorder', 'superfamily', 'family', 'genus', 'subgenus', 'group', 'subgroup',
        'Species', 'assembly_name', 'Sequence_ID', 'Contig', 'Start', 'End', 'Strand', 'Length',
        'motif_5prime_type', 'motif_3prime_type', 'longest_polyT_stretch', 'trailing_T_count',
        'noe_evalue', 'pol3_score', 'Sequence'
    ]

    with open(output_file, 'w', newline='') as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for sort_key, row in all_rows:
            writer.writerow(row)

    print(f"Exported {total_sequences} sequences from {species_count} species to {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Export all results from results/ folder to TSV',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Export with default settings
  python3 export_results_to_tsv.py -o results_summary.tsv

  # Specify custom paths
  python3 export_results_to_tsv.py --results-dir results --taxonomy data/genomes_annotated.tsv -o output.tsv
        """
    )

    parser.add_argument('--results-dir', default='results',
                       help='Path to results directory (default: results)')
    parser.add_argument('--taxonomy', default='data/genomes_annotated.tsv',
                       help='Path to taxonomy TSV file (default: data/genomes_annotated.tsv)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output TSV file')

    args = parser.parse_args()

    # Load taxonomy
    print(f"Loading taxonomy from {args.taxonomy}...", file=sys.stderr)
    taxonomy = load_taxonomy(args.taxonomy)
    print(f"Loaded taxonomy for {len(taxonomy)} species", file=sys.stderr)

    # Export results
    print(f"Exporting results from {args.results_dir}...", file=sys.stderr)
    export_results(args.results_dir, taxonomy, args.output)


if __name__ == "__main__":
    main()
