#!/usr/bin/env python3
"""
Script to concatenate all lncrna.pdf files in results/ folder
ordered by taxonomic lineage using pdftk
"""

import os
import sys
import subprocess
from pathlib import Path
import argparse

def get_taxonomic_lineage(species_name, tsv_file):
    """
    Extract taxonomic lineage for a given species from the TSV file.
    Returns a list of taxonomic levels for sorting.
    """
    try:
        with open(tsv_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 14 and fields[0] == species_name:
                    # Extract taxonomic levels for sorting
                    taxonomic_levels = []

                    # Order (index 5), Suborder, Infraorder, Superfamily, Family, Genus, Subgenus, Group, Subgroup
                    level_indices = [5, 6, 7, 8, 9, 10, 11, 12, 13]

                    for idx in level_indices:
                        if idx < len(fields):
                            level = fields[idx] if fields[idx] != 'NA' else ''
                            taxonomic_levels.append(level)
                        else:
                            taxonomic_levels.append('')

                    # Add species name at the end
                    taxonomic_levels.append(species_name)

                    return taxonomic_levels

        return None

    except Exception as e:
        print(f"Error reading taxonomy for {species_name}: {str(e)}", file=sys.stderr)
        return None

def collect_and_sort_pdfs(results_dir='results', tsv_file='data/genomes_annotated.tsv'):
    """
    Collect all lncrna.pdf files and sort them by taxonomic lineage.
    Returns a list of (pdf_path, taxonomic_lineage) tuples sorted by taxonomy.
    """
    pdf_files = []

    # Find all lncrna.pdf files
    results_path = Path(results_dir)
    for pdf_path in results_path.glob('*/lncrna.pdf'):
        # Extract species name from directory
        species_name = pdf_path.parent.name.replace('_', ' ')

        # Get taxonomic lineage
        lineage = get_taxonomic_lineage(species_name, tsv_file)

        if lineage:
            pdf_files.append((str(pdf_path), lineage, species_name))
        else:
            # If no taxonomy found, add with species name only for sorting
            print(f"Warning: No taxonomic lineage found for {species_name}", file=sys.stderr)
            pdf_files.append((str(pdf_path), [''] * 9 + [species_name], species_name))

    # Sort by taxonomic lineage
    pdf_files.sort(key=lambda x: x[1])

    return pdf_files

def concatenate_pdfs(pdf_files, output_file='all_lncrna.pdf'):
    """
    Use pdftk to concatenate PDFs in the specified order.
    """
    if not pdf_files:
        print("No PDF files found to concatenate.", file=sys.stderr)
        return False

    # Extract just the file paths
    pdf_paths = [pdf[0] for pdf in pdf_files]

    # Build pdftk command
    cmd = ['pdftk'] + pdf_paths + ['cat', 'output', output_file]

    print(f"Concatenating {len(pdf_paths)} PDF files...")
    print("Order of concatenation:")
    for i, (path, lineage, species) in enumerate(pdf_files, 1):
        lineage_str = ' > '.join([l for l in lineage[:-1] if l])
        if lineage_str:
            print(f"  {i:3d}. {species} ({lineage_str})")
        else:
            print(f"  {i:3d}. {species}")

    try:
        # Run pdftk
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"\nSuccessfully created {output_file}")

            # Get file size
            output_path = Path(output_file)
            if output_path.exists():
                size_mb = output_path.stat().st_size / (1024 * 1024)
                print(f"Output file size: {size_mb:.2f} MB")

            return True
        else:
            print(f"Error running pdftk: {result.stderr}", file=sys.stderr)
            return False

    except FileNotFoundError:
        print("Error: pdftk not found. Please install pdftk first.", file=sys.stderr)
        print("On Ubuntu/Debian: sudo apt-get install pdftk", file=sys.stderr)
        print("On macOS: brew install pdftk-java", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error concatenating PDFs: {str(e)}", file=sys.stderr)
        return False

def main():
    parser = argparse.ArgumentParser(
        description='Concatenate all lncrna.pdf files ordered by taxonomic lineage'
    )
    parser.add_argument(
        '-r', '--results-dir',
        default='results',
        help='Directory containing species subdirectories with lncrna.pdf files (default: results)'
    )
    parser.add_argument(
        '-t', '--tsv-file',
        default='data/genomes_annotated.tsv',
        help='Path to genomes_annotated.tsv file (default: data/genomes_annotated.tsv)'
    )
    parser.add_argument(
        '-o', '--output',
        default='all_lncrna.pdf',
        help='Output PDF filename (default: all_lncrna.pdf)'
    )

    args = parser.parse_args()

    # Check if results directory exists
    if not os.path.exists(args.results_dir):
        print(f"Error: Results directory '{args.results_dir}' not found.", file=sys.stderr)
        sys.exit(1)

    # Check if TSV file exists
    if not os.path.exists(args.tsv_file):
        print(f"Error: TSV file '{args.tsv_file}' not found.", file=sys.stderr)
        sys.exit(1)

    # Collect and sort PDFs
    pdf_files = collect_and_sort_pdfs(args.results_dir, args.tsv_file)

    if not pdf_files:
        print("No lncrna.pdf files found in the results directory.", file=sys.stderr)
        sys.exit(1)

    # Concatenate PDFs
    if concatenate_pdfs(pdf_files, args.output):
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()