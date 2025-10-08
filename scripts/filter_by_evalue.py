#!/usr/bin/env python3
"""
Filter TSV file by noe_evalue threshold and extract Sequence_IDs.

Usage:
    python3 filter_by_evalue.py results_summary_new.tsv
    python3 filter_by_evalue.py results_summary_new.tsv -e 1e-50 -o filtered_ids.txt
"""

import sys
import csv
import argparse
from pathlib import Path


def parse_evalue(evalue_str):
    """Parse e-value string, handling 'NA' and various formats"""
    if not evalue_str or evalue_str == 'NA' or evalue_str.strip() == '':
        return None

    try:
        return float(evalue_str)
    except ValueError:
        return None


def filter_by_evalue(tsv_file, threshold, output_file=None):
    """Filter TSV by e-value threshold and output Sequence_IDs"""

    if not Path(tsv_file).exists():
        print(f"Error: File not found: {tsv_file}", file=sys.stderr)
        sys.exit(1)

    # Determine output destination
    if output_file:
        try:
            outfile = open(output_file, 'w')
        except IOError as e:
            print(f"Error opening output file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        outfile = sys.stdout

    matching_ids = []
    total_count = 0
    filtered_count = 0
    na_count = 0

    try:
        with open(tsv_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')

            # Check if required columns exist
            if 'Sequence_ID' not in reader.fieldnames:
                print("Error: 'Sequence_ID' column not found in TSV", file=sys.stderr)
                sys.exit(1)
            if 'noe_evalue' not in reader.fieldnames:
                print("Error: 'noe_evalue' column not found in TSV", file=sys.stderr)
                sys.exit(1)

            for row in reader:
                total_count += 1
                seq_id = row['Sequence_ID']
                evalue_str = row['noe_evalue']

                evalue = parse_evalue(evalue_str)

                if evalue is None:
                    na_count += 1
                    continue

                if evalue < threshold:
                    matching_ids.append(seq_id)
                    outfile.write(f"{seq_id}\n")
                    filtered_count += 1

        # Print statistics to stderr
        print(f"Total sequences: {total_count}", file=sys.stderr)
        print(f"E-value < {threshold}: {filtered_count}", file=sys.stderr)
        print(f"NA or invalid e-values: {na_count}", file=sys.stderr)

    finally:
        if output_file:
            outfile.close()

    return matching_ids


def main():
    parser = argparse.ArgumentParser(
        description='Filter TSV file by noe_evalue threshold',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Filter with default threshold (1e-50)
  python3 filter_by_evalue.py results_summary_new.tsv

  # Filter with custom threshold
  python3 filter_by_evalue.py results_summary_new.tsv -e 1e-100

  # Save to output file
  python3 filter_by_evalue.py results_summary_new.tsv -e 1e-50 -o filtered_ids.txt
        """
    )

    parser.add_argument('input',
                       help='Input TSV file')
    parser.add_argument('-e', '--evalue', type=float, default=1e-50,
                       help='E-value threshold (default: 1e-50)')
    parser.add_argument('-o', '--output',
                       help='Output file for Sequence_IDs (default: stdout)')

    args = parser.parse_args()

    filter_by_evalue(args.input, args.evalue, args.output)


if __name__ == "__main__":
    main()
