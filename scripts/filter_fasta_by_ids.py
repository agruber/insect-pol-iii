#!/usr/bin/env python3
"""
Filter FASTA sequences by ID list.
Takes FASTA from stdin and ID list file as argument.

Usage:
    cat sequences.fa | python3 filter_fasta_by_ids.py ids.txt
    python3 filter_fasta_by_ids.py ids.txt < sequences.fa
"""

import sys
import argparse
from typing import Set


def load_ids(id_file: str) -> Set[str]:
    """Load sequence IDs from file (one per line)"""
    ids = set()

    try:
        with open(id_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    ids.add(line)
    except IOError as e:
        print(f"Error reading ID file: {e}", file=sys.stderr)
        sys.exit(1)

    return ids


def extract_id_from_header(header: str) -> str:
    """
    Extract sequence ID from FASTA header.
    Handles various formats:
    - >ID
    - >ID|other|info
    - >ID description
    """
    # Remove '>' if present
    if header.startswith('>'):
        header = header[1:]

    # Split by whitespace or pipe and take first part
    if ' ' in header:
        id_part = header.split()[0]
    else:
        id_part = header

    # If there's a pipe, take the first part before the pipe
    if '|' in id_part:
        id_part = id_part.split('|')[0]

    return id_part


def filter_fasta(input_stream, id_set: Set[str], output_stream):
    """Filter FASTA sequences based on ID set"""

    current_header = None
    current_seq = []
    matched_count = 0
    total_count = 0

    for line in input_stream:
        line = line.strip()

        if line.startswith('>'):
            # Process previous sequence if any
            if current_header:
                total_count += 1
                seq_id = extract_id_from_header(current_header)

                if seq_id in id_set:
                    # Write unwrapped FASTA
                    output_stream.write(f"{current_header}\n")
                    output_stream.write(f"{''.join(current_seq)}\n")
                    matched_count += 1

            # Start new sequence
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)

    # Don't forget the last sequence
    if current_header:
        total_count += 1
        seq_id = extract_id_from_header(current_header)

        if seq_id in id_set:
            output_stream.write(f"{current_header}\n")
            output_stream.write(f"{''.join(current_seq)}\n")
            matched_count += 1

    # Print statistics to stderr
    print(f"Total sequences in input: {total_count}", file=sys.stderr)
    print(f"Sequences matching IDs: {matched_count}", file=sys.stderr)
    print(f"IDs in list: {len(id_set)}", file=sys.stderr)

    if matched_count < len(id_set):
        print(f"Warning: {len(id_set) - matched_count} IDs from list not found in FASTA",
              file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Filter FASTA sequences by ID list',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Filter from stdin
  cat sequences.fa | python3 filter_fasta_by_ids.py ids.txt

  # Using input redirection
  python3 filter_fasta_by_ids.py ids.txt < sequences.fa

  # Chained with other commands
  cat sequences.fa | python3 filter_fasta_by_ids.py ids.txt > filtered.fa

Note: Outputs unwrapped FASTA (sequence on single line)
        """
    )

    parser.add_argument('ids',
                       help='File containing sequence IDs (one per line)')

    args = parser.parse_args()

    # Load ID list
    id_set = load_ids(args.ids)

    if not id_set:
        print("Error: No IDs loaded from file", file=sys.stderr)
        sys.exit(1)

    # Filter FASTA from stdin
    filter_fasta(sys.stdin, id_set, sys.stdout)


if __name__ == "__main__":
    main()
