#!/usr/bin/env python3
"""
Extract uppercase nucleotides, build PWMs for Pol II and Pol III,
and score sequences using MATCH algorithm.
Reads from STDIN and writes TSV scores to STDOUT.
"""

import sys
import yaml
import argparse
import math
import json
import tempfile
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

# Import MATCH algorithm
sys.path.insert(0, str(Path(__file__).parent))
from match_algorithm import MATCHAlgorithm


def load_config(config_file):
    """Load RNA polymerase classifications from config file."""
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        pol2_types = set(config['rna_polymerases']['pol2']['types'])
        pol3_types = set(config['rna_polymerases']['pol3']['types'])
        return pol2_types, pol3_types
    except Exception as e:
        print(f"Error loading config file: {e}", file=sys.stderr)
        sys.exit(1)


def extract_uppercase_nucleotides(sequence):
    """Extract only uppercase nucleotides from a sequence."""
    return ''.join([nt for nt in str(sequence) if nt.isupper()])


def classify_sequence(record_id, pol2_types, pol3_types):
    """Classify sequence as pol2, pol3, or unknown based on gene type."""
    # Extract gene type from record ID
    if '|' in record_id:
        gene_type = record_id.split('|')[1]
    elif '-' in record_id:
        parts = record_id.split('-')
        if len(parts) > 1:
            gene_type = parts[1]
        else:
            gene_type = "unknown"
    else:
        gene_type = "unknown"

    if gene_type in pol2_types:
        return "pol2"
    elif gene_type in pol3_types:
        return "pol3"
    else:
        return "unknown"


def build_pwm_from_sequences(sequences):
    """Build a Position Weight Matrix from a list of sequences."""
    if not sequences:
        return None

    # Check if all sequences have the same length
    seq_length = len(sequences[0])
    if not all(len(seq) == seq_length for seq in sequences):
        print(f"Warning: sequences have different lengths", file=sys.stderr)
        return None

    # Count nucleotides at each position
    counts = defaultdict(lambda: defaultdict(int))
    for seq in sequences:
        for i, nt in enumerate(seq.upper()):
            if nt in 'ACGT':
                counts[i][nt] += 1

    # Calculate frequencies and information content
    frequencies = {'A': [], 'C': [], 'G': [], 'T': []}
    information = []

    num_seqs = len(sequences)
    for i in range(seq_length):
        total = sum(counts[i].values())
        if total == 0:
            # No valid nucleotides at this position
            for nt in 'ACGT':
                frequencies[nt].append(0.25)  # Uniform background
            information.append(0.0)
        else:
            # Calculate frequencies
            freqs = {}
            for nt in 'ACGT':
                freqs[nt] = counts[i][nt] / total
                frequencies[nt].append(freqs[nt])

            # Calculate information content (bits)
            ic = 0
            for nt in 'ACGT':
                if freqs[nt] > 0:
                    ic += freqs[nt] * math.log2(freqs[nt] / 0.25)  # 0.25 = background freq
            information.append(max(0, ic))  # Ensure non-negative

    # Define core region (middle third of the motif)
    core_start = seq_length // 3
    core_end = seq_length - seq_length // 3

    pwm = {
        'length': seq_length,
        'frequencies': frequencies,
        'information': information,
        'core_start': core_start,
        'core_end': core_end,
        'num_sequences': num_seqs
    }

    return pwm


def build_split_pwm(all_sequences, pol2_sequences, pol3_sequences, split_position):
    """Build split PWM: first X positions from all sequences, rest separately by polymerase."""
    if not all_sequences or not pol2_sequences or not pol3_sequences:
        return None, None

    # Check sequence lengths
    seq_length = len(all_sequences[0])
    if split_position >= seq_length:
        print(f"Error: split position {split_position} >= sequence length {seq_length}", file=sys.stderr)
        return None, None

    # Build common prefix PWM from all sequences
    prefix_sequences = [seq[:split_position] for seq in all_sequences]
    common_prefix_pwm = build_pwm_from_sequences(prefix_sequences)

    if not common_prefix_pwm:
        return None, None

    # Build suffix PWMs separately for each polymerase
    pol2_suffix_sequences = [seq[split_position:] for seq in pol2_sequences]
    pol3_suffix_sequences = [seq[split_position:] for seq in pol3_sequences]

    pol2_suffix_pwm = build_pwm_from_sequences(pol2_suffix_sequences)
    pol3_suffix_pwm = build_pwm_from_sequences(pol3_suffix_sequences)

    if not pol2_suffix_pwm or not pol3_suffix_pwm:
        return None, None

    # Combine prefix + suffix for each polymerase
    def combine_pwms(prefix_pwm, suffix_pwm):
        combined_frequencies = {'A': [], 'C': [], 'G': [], 'T': []}
        combined_information = []

        # Add prefix
        for nt in 'ACGT':
            combined_frequencies[nt].extend(prefix_pwm['frequencies'][nt])
        combined_information.extend(prefix_pwm['information'])

        # Add suffix
        for nt in 'ACGT':
            combined_frequencies[nt].extend(suffix_pwm['frequencies'][nt])
        combined_information.extend(suffix_pwm['information'])

        total_length = prefix_pwm['length'] + suffix_pwm['length']

        return {
            'length': total_length,
            'frequencies': combined_frequencies,
            'information': combined_information,
            'core_start': total_length // 3,
            'core_end': total_length - total_length // 3,
            'num_sequences': suffix_pwm['num_sequences'],  # Use suffix count for polymerase-specific
            'split_position': split_position
        }

    pol2_pwm = combine_pwms(common_prefix_pwm, pol2_suffix_pwm)
    pol3_pwm = combine_pwms(common_prefix_pwm, pol3_suffix_pwm)

    return pol2_pwm, pol3_pwm


def main():
    parser = argparse.ArgumentParser(description='Build PWMs and score sequences using MATCH algorithm')
    parser.add_argument('-c', '--config', default='config.yaml',
                       help='Configuration file (default: config.yaml)')
    parser.add_argument('-e', '--exclude', action='append',
                       help='Exclude ncRNA types from PWM building (can be used multiple times, e.g., -e U1 -e U6)')
    parser.add_argument('--split', type=int,
                       help='Split mode: use first X nucleotides from all sequences, rest separately by polymerase')

    args = parser.parse_args()

    # Load RNA polymerase classifications
    pol2_types, pol3_types = load_config(args.config)

    # Set up exclusion list
    exclude_types = set(args.exclude) if args.exclude else set()

    if exclude_types:
        print(f"Excluding the following ncRNA types from PWM building: {', '.join(sorted(exclude_types))}", file=sys.stderr)

    # Group sequences by polymerase type
    pol2_sequences = []
    pol3_sequences = []
    all_sequences = []

    # Read FASTA from stdin
    for record in SeqIO.parse(sys.stdin, "fasta"):
        # Extract uppercase nucleotides
        uppercase_seq = extract_uppercase_nucleotides(record.seq)

        # Skip sequences without uppercase nucleotides
        if not uppercase_seq:
            continue

        # Classify sequence
        pol_type = classify_sequence(record.id, pol2_types, pol3_types)

        # Extract gene type for exclusion checking
        gene_type = record.id.split('|')[1] if '|' in record.id else record.id.split('-')[1] if '-' in record.id else "unknown"

        # Store all sequences for scoring
        all_sequences.append((record.id, uppercase_seq, pol_type))

        # Group sequences for PWM building (apply exclusions)
        if pol_type == "pol2" and gene_type not in exclude_types:
            pol2_sequences.append(uppercase_seq)
        elif pol_type == "pol3" and gene_type not in exclude_types:
            pol3_sequences.append(uppercase_seq)

    # Build PWMs
    if args.split:
        # Split mode: first X nucleotides from all sequences, rest separately
        all_pwm_sequences = pol2_sequences + pol3_sequences
        print(f"Split mode: Using first {args.split} nucleotides from all {len(all_pwm_sequences)} sequences...", file=sys.stderr)
        print(f"  Pol II suffix from {len(pol2_sequences)} sequences, Pol III suffix from {len(pol3_sequences)} sequences", file=sys.stderr)

        pol2_pwm, pol3_pwm = build_split_pwm(all_pwm_sequences, pol2_sequences, pol3_sequences, args.split)
    else:
        # Standard mode: separate PWMs
        print(f"Building Pol II PWM from {len(pol2_sequences)} sequences...", file=sys.stderr)
        pol2_pwm = build_pwm_from_sequences(pol2_sequences) if len(pol2_sequences) >= 2 else None

        print(f"Building Pol III PWM from {len(pol3_sequences)} sequences...", file=sys.stderr)
        pol3_pwm = build_pwm_from_sequences(pol3_sequences) if len(pol3_sequences) >= 2 else None

    if not pol2_pwm or not pol3_pwm:
        print("Error: Unable to build PWMs. Need at least 2 sequences per group.", file=sys.stderr)
        sys.exit(1)

    # Initialize MATCH algorithms
    pol2_matcher = MATCHAlgorithm(pol2_pwm)
    pol3_matcher = MATCHAlgorithm(pol3_pwm)

    # Output TSV header
    print("Sequence_ID\tPol2_Score\tPol3_Score")

    # Score all sequences
    for seq_id, seq, pol_type in all_sequences:
        # Calculate normalized MSS scores (0-1 range)
        try:
            pol2_score = pol2_matcher.calculate_mss(seq) if len(seq) >= pol2_matcher.length else 0.0
        except:
            pol2_score = 0.0

        try:
            pol3_score = pol3_matcher.calculate_mss(seq) if len(seq) >= pol3_matcher.length else 0.0
        except:
            pol3_score = 0.0

        print(f"{seq_id}\t{pol2_score:.2f}\t{pol3_score:.2f}")


if __name__ == '__main__':
    main()