#!/usr/bin/env python3
"""
Extract 5' and 3' regions from RNA sequences based on processed structure.
Only processes sequences with 'x' regions in their structure.

Usage:
    python3 extract_structure_regions.py -f sequences.fa -s structures.txt --5prime 2 --3prime 2
"""

import sys
import argparse
from pathlib import Path
from typing import Dict, Tuple, Optional


def parse_fasta(fasta_file: str) -> Dict[str, str]:
    """Parse FASTA file and return dict of header -> sequence"""
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                # Start new sequence
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget last sequence
        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences


def parse_structure_file(structure_file: str) -> Dict[str, Tuple[str, str]]:
    """
    Parse structure file and return dict of header -> (unprocessed_structure, processed_structure)
    """
    structures = {}

    with open(structure_file, 'r') as f:
        lines = [line.strip() for line in f]

    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            header = lines[i]
            if i + 2 < len(lines):
                unprocessed = lines[i + 1]
                processed = lines[i + 2]
                structures[header] = (unprocessed, processed)
                i += 3
            else:
                i += 1
        else:
            i += 1

    return structures


def find_x_region(structure: str) -> Optional[Tuple[int, int]]:
    """
    Find the x region boundaries (start and end indices).
    Returns None if no 'x' found.
    """
    if 'x' not in structure:
        return None

    # Find first and last 'x'
    first_x = structure.index('x')
    last_x = structure.rindex('x')

    return (first_x, last_x)


def extract_regions(sequence: str, structure: str, extend_5prime: int, extend_3prime: int) -> Optional[str]:
    """
    Extract 5' and 3' regions with extensions.
    Returns the joined sequence with NNN separator, or None if no 'x' region.
    """
    x_region = find_x_region(structure)

    if not x_region:
        return None

    x_start, x_end = x_region

    # 5' region: from start to first 'x', then extend INTO x region
    fivep_start = 0
    fivep_end = x_start + extend_5prime  # Extend into x region

    # 3' region: from last 'x' to end, but start earlier (extend INTO x region)
    threep_start = x_end + 1 - extend_3prime  # Start earlier by extending into x region
    threep_end = len(sequence)

    # Ensure indices are within bounds
    fivep_end = min(fivep_end, len(sequence))
    threep_start = max(0, threep_start)

    # Extract regions
    fivep_region = sequence[fivep_start:fivep_end]
    threep_region = sequence[threep_start:threep_end]

    # Join with NNN
    result = fivep_region + 'NNN' + threep_region

    return result


def main():
    parser = argparse.ArgumentParser(
        description='Extract 5\' and 3\' regions from sequences with structured cores',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract with 2nt extensions on both ends
  python3 extract_structure_regions.py -f lncrna.fa -s lncrna.fa.structure --5prime 2 --3prime 2

  # No extension
  python3 extract_structure_regions.py -f lncrna.fa -s lncrna.fa.structure --5prime 0 --3prime 0

  # Different extensions
  python3 extract_structure_regions.py -f lncrna.fa -s lncrna.fa.structure --5prime 3 --3prime 1 -o output.fa
        """
    )

    parser.add_argument('-f', '--fasta', required=True,
                       help='Input FASTA file')
    parser.add_argument('-s', '--structure', required=True,
                       help='Structure file (output from analyze_structure.py)')
    parser.add_argument('--5prime', type=int, default=0, dest='extend_5prime',
                       help='Number of nucleotides to extend 5\' region into x region (default: 0)')
    parser.add_argument('--3prime', type=int, default=0, dest='extend_3prime',
                       help='Number of nucleotides to extend 3\' region into x region (default: 0)')
    parser.add_argument('-o', '--output',
                       help='Output FASTA file (default: stdout)')

    args = parser.parse_args()

    # Check input files exist
    if not Path(args.fasta).exists():
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    if not Path(args.structure).exists():
        print(f"Error: Structure file not found: {args.structure}", file=sys.stderr)
        sys.exit(1)

    # Parse input files
    sequences = parse_fasta(args.fasta)
    structures = parse_structure_file(args.structure)

    # Determine output destination
    if args.output:
        try:
            outfile = open(args.output, 'w')
        except IOError as e:
            print(f"Error opening output file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        outfile = sys.stdout

    # Process sequences
    processed_count = 0
    skipped_count = 0

    try:
        for header in sequences:
            if header not in structures:
                print(f"Warning: No structure found for {header}", file=sys.stderr)
                skipped_count += 1
                continue

            sequence = sequences[header]
            unprocessed_struct, processed_struct = structures[header]

            # Extract regions
            result_seq = extract_regions(sequence, processed_struct,
                                        args.extend_5prime, args.extend_3prime)

            if result_seq:
                # Write without line wrapping
                outfile.write(f"{header}\n")
                outfile.write(f"{result_seq}\n")
                processed_count += 1
            else:
                skipped_count += 1

        print(f"Processed {processed_count} sequences, skipped {skipped_count} (no 'x' region)",
              file=sys.stderr)

    finally:
        if args.output:
            outfile.close()


if __name__ == "__main__":
    main()
