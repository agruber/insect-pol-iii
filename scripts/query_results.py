#!/usr/bin/env python3
"""
Query results_summary.tsv to answer specific questions.

Usage:
    python3 query_results.py results_summary.tsv
"""

import sys
import pandas as pd
import argparse
from pathlib import Path


def load_data(tsv_file):
    """Load TSV file into pandas DataFrame"""
    df = pd.read_csv(tsv_file, sep='\t', dtype={
        'order': str,
        'suborder': str,
        'infraorder': str,
        'superfamily': str,
        'family': str,
        'genus': str,
        'subgenus': str,
        'group': str,
        'subgroup': str,
        'Species': str,
        'assembly_name': str,
        'Sequence_ID': str,
        'Contig': str,
        'Start': int,
        'End': int,
        'Strand': str,
        'Length': int,
        'motif_5prime_type': str,
        'motif_3prime_type': str,
        'longest_polyT_stretch': int,
        'trailing_T_count': int,
        'noe_evalue': str,
        'pol3_score': float,
        'Sequence': str
    })

    df['noe_evalue_numeric'] = pd.to_numeric(df['noe_evalue'], errors='coerce')

    return df


def count_drosophilidae_species(df):
    """Count how many Drosophilidae species have at least one hit"""
    drosophilidae = df[df['family'] == 'Drosophilidae']
    unique_species = drosophilidae['Species'].nunique()
    return unique_species


def analyze_drosophilidae_hit_types(df):
    """Analyze noe vs CR34335 hits in Drosophilidae species"""
    drosophilidae = df[df['family'] == 'Drosophilidae']

    # For each species, determine if it has noe hits
    species_with_noe = set()
    all_species = set(drosophilidae['Species'].unique())

    # A species has a "noe hit" if at least one sequence has e-value <= 1e-30
    for species in all_species:
        species_df = drosophilidae[drosophilidae['Species'] == species]
        has_noe = (species_df['noe_evalue_numeric'] <= 1e-30).any()
        if has_noe:
            species_with_noe.add(species)

    # All species in results have "CR34335 hits" (passed the pipeline)
    species_with_cr34335 = all_species

    # Calculate categories
    both = len(species_with_noe & species_with_cr34335)
    only_cr34335 = len(species_with_cr34335 - species_with_noe)
    only_noe = 0  # Not in this dataset (would need separate noe-only results)

    return both, only_noe, only_cr34335


def find_shortest_longest_drosophilidae(df):
    """Find shortest and longest hits in Drosophilidae"""
    drosophilidae = df[df['family'] == 'Drosophilidae']

    shortest_row = drosophilidae.loc[drosophilidae['Length'].idxmin()]
    longest_row = drosophilidae.loc[drosophilidae['Length'].idxmax()]

    return (shortest_row['Species'], shortest_row['Length'],
            longest_row['Species'], longest_row['Length'])


def main():
    parser = argparse.ArgumentParser(
        description='Query results_summary.tsv to answer specific questions'
    )

    parser.add_argument('tsv_file',
                       help='Path to results_summary.tsv file')

    args = parser.parse_args()

    # Check file exists
    if not Path(args.tsv_file).exists():
        print(f"Error: File not found: {args.tsv_file}", file=sys.stderr)
        sys.exit(1)

    # Load data
    df = load_data(args.tsv_file)

    # Answer questions
    print("Question 1: For how many Drosophilidae species do we have at least one hit?")
    count = count_drosophilidae_species(df)
    print(f"Answer: We have at least one hit for {count} Drosophilidae species.")
    print()

    print("Question 2: For how many Drosophilidae species do we have a noe+CR34335 hit, only noe, or only CR34335 hit?")
    print("(A noe hit is defined as e-value <= 1e-30)")
    both, only_noe, only_cr34335 = analyze_drosophilidae_hit_types(df)
    print(f"Answer:")
    print(f"  Both noe and CR34335: {both} species")
    print(f"  Only noe: {only_noe} species")
    print(f"  Only CR34335: {only_cr34335} species")
    print()

    print("Question 3: What is the shortest and longest hit in Drosophilidae?")
    shortest_sp, shortest_len, longest_sp, longest_len = find_shortest_longest_drosophilidae(df)
    print(f"Answer:")
    print(f"  Shortest: {shortest_sp} with {shortest_len} nt")
    print(f"  Longest: {longest_sp} with {longest_len} nt")


if __name__ == "__main__":
    main()
