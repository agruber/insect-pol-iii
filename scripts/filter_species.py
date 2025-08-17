#!/usr/bin/env python3
"""
Species filter script - outputs filtered species data
Filters genomes_annotated.tsv based on taxonomic queries and outputs selected columns
"""

import sys
import pandas as pd
import argparse

def clean_species_name(species_name):
    """Convert species name to filesystem-safe format"""
    return species_name.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '_')

def main():
    parser = argparse.ArgumentParser(description='Filter species based on taxonomic queries')
    parser.add_argument('tsv_file', help='Path to the genomes TSV file')
    parser.add_argument('query', help='Taxonomic query (species name or taxonomic group)')
    parser.add_argument('-t', '--taxonomy-level', choices=[
        'phylum', 'class', 'order', 'suborder', 'infraorder', 
        'superfamily', 'family', 'genus', 'subgenus', 
        'species_group', 'species_subgroup'
    ], help='Specify taxonomy level for the query')
    parser.add_argument('-f', '--format', choices=['tsv', 'names'], default='tsv',
                       help='Output format: tsv (full data) or names (cleaned species names only)')
    
    args = parser.parse_args()
    
    # Column mapping
    col_mapping = {
        'species': 0,
        'assembly_type': 1,
        'url': 2,
        'phylum': 3,
        'class': 4,
        'order': 5,
        'suborder': 6,
        'infraorder': 7,
        'superfamily': 8,
        'family': 9,
        'genus': 10,
        'subgenus': 11,
        'species_group': 12,
        'species_subgroup': 13
    }
    
    # Read the TSV file
    df = pd.read_csv(args.tsv_file, sep='\t', header=None)
    
    # Filter based on query
    if args.taxonomy_level:
        col_idx = col_mapping[args.taxonomy_level]
        filtered_df = df[df.iloc[:, col_idx].str.contains(args.query, case=False, na=False)]
    else:
        # Search in species name (first column) by default
        filtered_df = df[df.iloc[:, 0].str.contains(args.query, case=False, na=False)]
    
    if filtered_df.empty:
        print(f"No matches found for '{args.query}'", file=sys.stderr)
        sys.exit(1)
    
    # Output based on format
    if args.format == 'names':
        # Output only cleaned species names
        for _, row in filtered_df.iterrows():
            print(clean_species_name(row.iloc[0]))
    else:
        # Output full TSV data
        filtered_df.to_csv(sys.stdout, sep='\t', header=False, index=False)

if __name__ == "__main__":
    main()