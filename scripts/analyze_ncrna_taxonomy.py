#!/usr/bin/env python3
"""
Analyze ncRNA distribution across taxonomic groups.
Identifies missing ncRNAs at different taxonomic levels.
"""

import sys
import os
from pathlib import Path
from collections import defaultdict
import pandas as pd
import argparse
from typing import Dict, Set, List
import yaml

def load_taxonomy_data(genomes_file: str = 'data/genomes_annotated.tsv') -> pd.DataFrame:
    """Load taxonomy data from the annotated genomes file"""
    try:
        # Read file with variable column counts
        df = pd.read_csv(genomes_file, sep='\t', header=None)

        # Expected columns (14 total based on file structure)
        column_names = ['species', 'assembly_level', 'url', 'phylum', 'class',
                       'order', 'suborder', 'infraorder', 'superfamily',
                       'family', 'genus', 'subgenus', 'extra1', 'extra2']

        # Only assign columns that exist
        num_cols = min(len(column_names), len(df.columns))
        df.columns = column_names[:num_cols] + [f'extra_{i}' for i in range(num_cols, len(df.columns))]

        # Fill in missing columns with NaN
        for col in column_names:
            if col not in df.columns:
                df[col] = pd.NA

        # Convert species names to match directory format (spaces to underscores)
        df['species_dir'] = df['species'].str.replace(' ', '_')
        return df
    except Exception as e:
        print(f"Error loading taxonomy data: {e}", file=sys.stderr)
        return pd.DataFrame()

def load_polymerase_config(config_file: str = 'config.yaml') -> Dict:
    """Load RNA polymerase classification from config"""
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        return config.get('rna_polymerases', {})
    except Exception as e:
        print(f"Error loading config: {e}", file=sys.stderr)
        return {}

def get_species_ncrnas(species_dir: Path) -> Set[str]:
    """Extract ncRNA types found for a species from upstream.fa"""
    ncrna_types = set()
    upstream_file = species_dir / 'upstream.fa'

    if not upstream_file.exists():
        return ncrna_types

    try:
        with open(upstream_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Format: >species-id|ncRNA_type|...
                    parts = line[1:].strip().split('|')
                    if len(parts) >= 2:
                        ncrna_types.add(parts[1])
    except Exception:
        pass

    return ncrna_types

def analyze_by_taxonomy_level(df: pd.DataFrame, genomes_dir: Path,
                             level: str, polymerase_config: Dict) -> pd.DataFrame:
    """Analyze ncRNA distribution at a specific taxonomic level"""

    # Get all possible ncRNA types from config
    all_pol2 = set(polymerase_config.get('pol2', {}).get('types', []))
    all_pol3 = set(polymerase_config.get('pol3', {}).get('types', []))
    all_ncrnas = all_pol2 | all_pol3

    # Group by taxonomic level
    results = []

    for group_name, group_df in df.groupby(level):
        if pd.isna(group_name) or group_name == 'NA':
            continue

        # Collect ncRNAs for all species in this group
        group_ncrnas = defaultdict(int)  # ncRNA -> count of species with it
        species_count = 0
        species_with_data = 0

        for _, row in group_df.iterrows():
            species_dir = genomes_dir / row['species_dir']
            if species_dir.exists():
                species_count += 1
                species_ncrnas = get_species_ncrnas(species_dir)
                if species_ncrnas:
                    species_with_data += 1
                    for ncrna in species_ncrnas:
                        group_ncrnas[ncrna] += 1

        if species_count == 0:
            continue

        # Find missing ncRNAs
        found_ncrnas = set(group_ncrnas.keys())
        missing_ncrnas = all_ncrnas - found_ncrnas
        missing_pol2 = all_pol2 - found_ncrnas
        missing_pol3 = all_pol3 - found_ncrnas

        # Calculate coverage
        coverage_pct = (len(found_ncrnas) / len(all_ncrnas) * 100) if all_ncrnas else 0

        results.append({
            level: group_name,
            'n_species': species_count,
            'n_with_data': species_with_data,
            'n_ncrnas_found': len(found_ncrnas),
            'n_ncrnas_total': len(all_ncrnas),
            'coverage_pct': coverage_pct,
            'n_missing': len(missing_ncrnas),
            'missing_pol2': ','.join(sorted(missing_pol2)) if missing_pol2 else '-',
            'missing_pol3': ','.join(sorted(missing_pol3)) if missing_pol3 else '-',
            'found_ncrnas': ','.join(sorted(found_ncrnas)) if found_ncrnas else '-'
        })

    return pd.DataFrame(results)

def analyze_rare_ncrnas(df: pd.DataFrame, genomes_dir: Path) -> pd.DataFrame:
    """Find ncRNAs that are rare across all species"""
    ncrna_counts = defaultdict(int)
    total_species_with_data = 0

    for _, row in df.iterrows():
        species_dir = genomes_dir / row['species_dir']
        if species_dir.exists():
            species_ncrnas = get_species_ncrnas(species_dir)
            if species_ncrnas:
                total_species_with_data += 1
                for ncrna in species_ncrnas:
                    ncrna_counts[ncrna] += 1

    # Create dataframe of ncRNA frequencies
    rare_data = []
    for ncrna, count in sorted(ncrna_counts.items(), key=lambda x: x[1]):
        percentage = (count / total_species_with_data * 100) if total_species_with_data > 0 else 0
        rare_data.append({
            'ncRNA': ncrna,
            'species_count': count,
            'percentage': percentage,
            'status': 'rare' if percentage < 10 else ('uncommon' if percentage < 50 else 'common')
        })

    return pd.DataFrame(rare_data)

def main():
    parser = argparse.ArgumentParser(
        description='Analyze ncRNA distribution across taxonomic groups'
    )
    parser.add_argument('-d', '--directory', default='genomes',
                       help='Genomes directory (default: genomes)')
    parser.add_argument('-t', '--taxonomy', default='data/genomes_annotated.tsv',
                       help='Taxonomy file (default: data/genomes_annotated.tsv)')
    parser.add_argument('-c', '--config', default='config.yaml',
                       help='Configuration file (default: config.yaml)')
    parser.add_argument('-l', '--level',
                       choices=['genus', 'family', 'order', 'class', 'all'],
                       default='all',
                       help='Taxonomic level to analyze (default: all)')
    parser.add_argument('-o', '--output', default=None,
                       help='Output file (default: stdout)')
    parser.add_argument('--rare', action='store_true',
                       help='Show rare ncRNAs analysis')

    args = parser.parse_args()

    genomes_dir = Path(args.directory)
    if not genomes_dir.exists():
        print(f"Error: Directory {genomes_dir} not found", file=sys.stderr)
        sys.exit(1)

    # Load data
    print("Loading taxonomy data...", file=sys.stderr)
    df = load_taxonomy_data(args.taxonomy)
    if df.empty:
        print("Error: Could not load taxonomy data", file=sys.stderr)
        sys.exit(1)

    polymerase_config = load_polymerase_config(args.config)

    output_parts = []

    # Analyze rare ncRNAs if requested
    if args.rare:
        print("Analyzing rare ncRNAs...", file=sys.stderr)
        rare_df = analyze_rare_ncrnas(df, genomes_dir)
        output_parts.append("=== Rare ncRNA Analysis ===\n")
        output_parts.append(rare_df.to_string(index=False))
        output_parts.append("\n\n")

    # Analyze by taxonomic level
    levels = ['genus', 'family', 'order', 'class'] if args.level == 'all' else [args.level]

    for level in levels:
        print(f"Analyzing at {level} level...", file=sys.stderr)
        level_df = analyze_by_taxonomy_level(df, genomes_dir, level, polymerase_config)

        if not level_df.empty:
            # Sort by coverage percentage
            level_df = level_df.sort_values('coverage_pct', ascending=True)

            output_parts.append(f"=== Analysis by {level.upper()} ===\n")

            # Show groups with missing ncRNAs
            missing_df = level_df[level_df['n_missing'] > 0]
            if not missing_df.empty:
                output_parts.append(f"\n{level.capitalize()} groups with missing ncRNAs:\n")
                output_parts.append(missing_df.to_string(index=False))
            else:
                output_parts.append(f"\nAll {level} groups have complete ncRNA coverage!")

            # Summary statistics
            output_parts.append(f"\n\nSummary for {level}:")
            output_parts.append(f"  Total groups: {len(level_df)}")
            output_parts.append(f"  Groups with complete coverage: {len(level_df[level_df['n_missing'] == 0])}")
            output_parts.append(f"  Average coverage: {level_df['coverage_pct'].mean():.1f}%")
            output_parts.append("\n\n")

    # Combine output
    output = ''.join(output_parts)

    # Write output
    if args.output:
        with open(args.output, 'w') as f:
            f.write(output)
        print(f"Output written to {args.output}", file=sys.stderr)
    else:
        print(output)

if __name__ == "__main__":
    main()