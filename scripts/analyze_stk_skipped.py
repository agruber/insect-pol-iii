#!/usr/bin/env python3
"""
Analyze STK file to determine how many nucleotides were skipped where NNN appears.

Usage:
    python3 analyze_stk_skipped.py structures/CR34335_melanogaster_subgroup.stk results_Oct_3_25/results_summary.tsv
"""

import sys
import pandas as pd
from pathlib import Path


def parse_stk_file(stk_file):
    """Parse Stockholm format file to extract sequences"""
    sequences = {}

    with open(stk_file, 'r') as f:
        for line in f:
            line = line.rstrip()

            # Skip comment lines, blank lines, and markup lines
            if not line or line.startswith('#') or line.startswith('//'):
                continue

            # Parse sequence line
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                seq = parts[1]

                # Append sequence (STK can have sequences on multiple lines)
                if seq_id in sequences:
                    sequences[seq_id] += seq
                else:
                    sequences[seq_id] = seq

    return sequences


def load_original_sequences(tsv_file):
    """Load original sequences from results_summary.tsv"""
    df = pd.read_csv(tsv_file, sep='\t')

    # Create a dictionary: Sequence_ID -> Sequence
    seq_dict = {}
    for _, row in df.iterrows():
        seq_id = row['Sequence_ID']
        seq = row['Sequence'].upper()  # Make uppercase for comparison
        seq_dict[seq_id] = seq

    return seq_dict


def find_skipped_nucleotides(stk_seq, original_seq):
    """
    Find how many nucleotides were skipped where NNN appears.

    stk_seq: sequence from STK file (with gaps and NNN)
    original_seq: original full sequence

    Returns: number of skipped nucleotides, or None if match failed
    """
    # Remove gaps from STK sequence
    stk_seq_nogaps = stk_seq.replace('-', '')

    # Convert U to T for comparison
    stk_seq_nogaps = stk_seq_nogaps.replace('U', 'T')

    # Find NNN position
    nnn_pos = stk_seq_nogaps.find('NNN')
    if nnn_pos == -1:
        return None

    # Extract 5' and 3' parts
    part_5prime = stk_seq_nogaps[:nnn_pos]
    part_3prime = stk_seq_nogaps[nnn_pos + 3:]  # Skip the NNN

    # Find these parts in the original sequence
    pos_5prime = original_seq.find(part_5prime)
    if pos_5prime == -1:
        return None

    # The 3' part should be after the 5' part
    search_start = pos_5prime + len(part_5prime)
    pos_3prime = original_seq.find(part_3prime, search_start)
    if pos_3prime == -1:
        return None

    # Calculate skipped nucleotides
    # Position where 5' ends
    end_5prime = pos_5prime + len(part_5prime)
    # Position where 3' starts
    start_3prime = pos_3prime

    # Number of nucleotides between them
    skipped = start_3prime - end_5prime

    return skipped


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 analyze_stk_skipped.py <stk_file> <results_summary.tsv>")
        sys.exit(1)

    stk_file = sys.argv[1]
    tsv_file = sys.argv[2]

    # Check files exist
    if not Path(stk_file).exists():
        print(f"Error: STK file not found: {stk_file}", file=sys.stderr)
        sys.exit(1)

    if not Path(tsv_file).exists():
        print(f"Error: TSV file not found: {tsv_file}", file=sys.stderr)
        sys.exit(1)

    # Parse STK file
    stk_sequences = parse_stk_file(stk_file)

    # Load original sequences
    original_sequences = load_original_sequences(tsv_file)

    # Analyze each sequence
    print("Sequence_ID\tSkipped_Nucleotides")

    skipped_values = []

    for seq_id, stk_seq in stk_sequences.items():
        # Try to find matching original sequence
        # The seq_id in STK might be shortened, try to find a match
        matched = False

        for orig_id, orig_seq in original_sequences.items():
            # Check if the STK seq_id is a prefix of the original ID
            # or if they match in some way
            if seq_id in orig_id or orig_id.startswith(seq_id.split('_')[0]):
                skipped = find_skipped_nucleotides(stk_seq, orig_seq)
                if skipped is not None:
                    print(f"{seq_id}\t{skipped}")
                    skipped_values.append(skipped)
                    matched = True
                    break

        if not matched:
            # Try direct match
            if seq_id in original_sequences:
                skipped = find_skipped_nucleotides(stk_seq, original_sequences[seq_id])
                if skipped is not None:
                    print(f"{seq_id}\t{skipped}")
                    skipped_values.append(skipped)
                else:
                    print(f"{seq_id}\tFailed to match", file=sys.stderr)
            else:
                print(f"{seq_id}\tNo original sequence found", file=sys.stderr)

    # Print statistics
    if skipped_values:
        print()
        print(f"Min: {min(skipped_values)} nt")
        print(f"Max: {max(skipped_values)} nt")


if __name__ == "__main__":
    main()
