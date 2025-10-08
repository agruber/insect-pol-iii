#!/usr/bin/env python3
"""
Analyze results_summary.tsv file with noeCR34335 lncRNA data.

Usage:
    python3 analyze_results_summary.py results_Oct_3_25/results_summary.tsv
"""

import sys
import pandas as pd
import argparse
from pathlib import Path
from collections import Counter


def load_data(tsv_file):
    """Load TSV file into pandas DataFrame"""
    print(f"Loading data from {tsv_file}...", file=sys.stderr)

    # Read TSV with proper handling of empty values
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
        'noe_evalue': str,  # Keep as string initially due to 'NA' values
        'pol3_score': float,
        'Sequence': str
    })

    # Convert noe_evalue to numeric, coercing 'NA' to NaN
    df['noe_evalue_numeric'] = pd.to_numeric(df['noe_evalue'], errors='coerce')

    print(f"Loaded {len(df)} sequences from {df['Species'].nunique()} species", file=sys.stderr)

    return df


def basic_statistics(df):
    """Print basic statistics about the dataset"""
    print("\n" + "="*80)
    print("BASIC STATISTICS")
    print("="*80)

    print(f"\nTotal sequences: {len(df)}")
    print(f"Total species: {df['Species'].nunique()}")
    print(f"Unique families: {df['family'].nunique()}")
    print(f"Unique superfamilies: {df['superfamily'].nunique()}")

    print(f"\nSequence length statistics:")
    print(f"  Min: {df['Length'].min()} nt")
    print(f"  Max: {df['Length'].max()} nt")
    print(f"  Mean: {df['Length'].mean():.1f} nt")
    print(f"  Median: {df['Length'].median():.1f} nt")


def motif_analysis(df):
    """Analyze motif combinations"""
    print("\n" + "="*80)
    print("MOTIF ANALYSIS")
    print("="*80)

    print("\n5' Motif distribution:")
    motif_5_counts = df['motif_5prime_type'].value_counts()
    for motif, count in motif_5_counts.items():
        pct = 100 * count / len(df)
        print(f"  {motif}: {count} ({pct:.1f}%)")

    print("\n3' Motif distribution:")
    motif_3_counts = df['motif_3prime_type'].value_counts()
    for motif, count in motif_3_counts.items():
        pct = 100 * count / len(df)
        print(f"  {motif}: {count} ({pct:.1f}%)")

    print("\nMotif pair combinations (top 10):")
    df['motif_pair'] = df['motif_5prime_type'] + ' + ' + df['motif_3prime_type']
    motif_pairs = df['motif_pair'].value_counts().head(10)
    for pair, count in motif_pairs.items():
        pct = 100 * count / len(df)
        print(f"  {pair}: {count} ({pct:.1f}%)")


def score_analysis(df):
    """Analyze scoring distributions"""
    print("\n" + "="*80)
    print("SCORE ANALYSIS")
    print("="*80)

    print("\nPol III promoter scores:")
    print(f"  Mean: {df['pol3_score'].mean():.3f}")
    print(f"  Median: {df['pol3_score'].median():.3f}")
    print(f"  Min: {df['pol3_score'].min():.3f}")
    print(f"  Max: {df['pol3_score'].max():.3f}")

    # Score bins
    print("\nPol III score distribution:")
    bins = [0, 0.8, 0.85, 0.9, 0.95, 1.0]
    labels = ['0.00-0.80', '0.80-0.85', '0.85-0.90', '0.90-0.95', '0.95-1.00']
    df['pol3_bin'] = pd.cut(df['pol3_score'], bins=bins, labels=labels, include_lowest=True)
    bin_counts = df['pol3_bin'].value_counts().sort_index()
    for bin_label, count in bin_counts.items():
        pct = 100 * count / len(df)
        print(f"  {bin_label}: {count} ({pct:.1f}%)")

    print("\nnoeCR34335 E-value statistics:")
    valid_evalues = df[df['noe_evalue_numeric'].notna()]
    na_count = df['noe_evalue_numeric'].isna().sum()

    print(f"  Valid E-values: {len(valid_evalues)}")
    print(f"  NA values: {na_count}")

    if len(valid_evalues) > 0:
        print(f"  Min E-value: {valid_evalues['noe_evalue_numeric'].min():.2e}")
        print(f"  Max E-value: {valid_evalues['noe_evalue_numeric'].max():.2e}")
        print(f"  Median E-value: {valid_evalues['noe_evalue_numeric'].median():.2e}")

        # Count sequences with very significant E-values
        very_sig = (valid_evalues['noe_evalue_numeric'] < 1e-50).sum()
        sig = ((valid_evalues['noe_evalue_numeric'] >= 1e-50) &
               (valid_evalues['noe_evalue_numeric'] < 1e-10)).sum()

        print(f"\n  E-value < 1e-50: {very_sig} ({100*very_sig/len(df):.1f}%)")
        print(f"  1e-50 ≤ E-value < 1e-10: {sig} ({100*sig/len(df):.1f}%)")


def taxonomy_distribution(df):
    """Analyze taxonomic distribution"""
    print("\n" + "="*80)
    print("TAXONOMIC DISTRIBUTION")
    print("="*80)

    print("\nTop 10 superfamilies by sequence count:")
    superfam_counts = df['superfamily'].value_counts().head(10)
    for superfam, count in superfam_counts.items():
        species_count = df[df['superfamily'] == superfam]['Species'].nunique()
        print(f"  {superfam}: {count} sequences from {species_count} species")

    print("\nTop 10 families by sequence count:")
    fam_counts = df['family'].value_counts().head(10)
    for fam, count in fam_counts.items():
        species_count = df[df['family'] == fam]['Species'].nunique()
        print(f"  {fam}: {count} sequences from {species_count} species")

    print("\nSequences per species (top 10):")
    species_counts = df['Species'].value_counts().head(10)
    for species, count in species_counts.items():
        print(f"  {species}: {count} sequences")


def polyt_analysis(df):
    """Analyze poly-T characteristics"""
    print("\n" + "="*80)
    print("POLY-T ANALYSIS")
    print("="*80)

    print("\nLongest poly-T stretch distribution:")
    polyt_counts = df['longest_polyT_stretch'].value_counts().sort_index()
    for stretch, count in polyt_counts.items():
        pct = 100 * count / len(df)
        print(f"  {stretch} T's: {count} ({pct:.1f}%)")

    print("\nTrailing T count statistics:")
    print(f"  Mean: {df['trailing_T_count'].mean():.1f}")
    print(f"  Median: {df['trailing_T_count'].median():.1f}")
    print(f"  Min: {df['trailing_T_count'].min()}")
    print(f"  Max: {df['trailing_T_count'].max()}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze noeCR34335 results summary TSV file',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('tsv_file',
                       help='Path to results_summary.tsv file')
    parser.add_argument('--export-filtered',
                       help='Export filtered data to file (e.g., high-confidence sequences)')

    args = parser.parse_args()

    # Check file exists
    if not Path(args.tsv_file).exists():
        print(f"Error: File not found: {args.tsv_file}", file=sys.stderr)
        sys.exit(1)

    # Load data
    df = load_data(args.tsv_file)

    # Run analyses
    basic_statistics(df)
    motif_analysis(df)
    score_analysis(df)
    taxonomy_distribution(df)
    polyt_analysis(df)

    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
