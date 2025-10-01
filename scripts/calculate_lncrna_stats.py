#!/usr/bin/env python3
"""
Calculate statistics for lncRNA sequences including motif presence and poly-T stretches.

Usage:
    python3 calculate_lncrna_stats.py input.fa > output.stats
    cat input.fa | python3 calculate_lncrna_stats.py > output.stats
"""

import sys
import re
import yaml
from Bio import SeqIO
from bioseq_lib import create_standard_parser, smart_open

def find_longest_polyt_between_motifs(sequence, start_pos, end_pos):
    """
    Find the longest stretch of consecutive Ts between two positions.

    Args:
        sequence: The sequence string (uppercase)
        start_pos: Start position (0-based, inclusive)
        end_pos: End position (0-based, exclusive)

    Returns:
        Length of the longest poly-T stretch
    """
    if start_pos >= end_pos or start_pos < 0 or end_pos > len(sequence):
        return 0

    # Extract the region between the motifs
    region = sequence[start_pos:end_pos]

    # Find all stretches of consecutive Ts
    poly_t_matches = re.finditer(r'T+', region)

    # Find the longest stretch
    max_length = 0
    for match in poly_t_matches:
        length = len(match.group())
        if length > max_length:
            max_length = length

    return max_length

def count_trailing_ts(sequence):
    """
    Count the number of T nucleotides at the 3' end.

    Args:
        sequence: The sequence string (uppercase)

    Returns:
        Number of trailing Ts
    """
    match = re.search(r'T+$', sequence)
    if match:
        return len(match.group())
    return 0

def trim_polyt_proper(sequence):
    """
    Properly trim polyT from 3' end, including cases like ...ATCGCATTTTTT.

    Steps:
    1. Trim trailing Ts
    2. If there's a single non-T at the end, check if there are Ts before it
    3. If yes, remove the non-T and continue trimming
    4. Repeat until sequence ends with multiple non-Ts or is all Ts

    Args:
        sequence: The sequence string (uppercase)

    Returns:
        Trimmed sequence
    """
    seq = sequence

    while len(seq) > 0:
        # Remove trailing Ts
        seq_no_t = seq.rstrip('T')

        if len(seq_no_t) == 0:
            # All Ts - return empty or the string based on preference
            return seq_no_t

        if len(seq_no_t) == len(seq):
            # No trailing Ts - we're done
            return seq

        # We removed some Ts. Now check if there's a single non-T at the end
        # If the last character is not T and there are Ts before it, remove it and continue
        if len(seq_no_t) >= 2 and seq_no_t[-1] != 'T' and seq_no_t[-2] == 'T':
            # Pattern like ...TN (where N is the last non-T)
            # Remove the last non-T and continue
            seq = seq_no_t[:-1]
        else:
            # Either ends with multiple non-Ts or just one non-T with no Ts before
            # We're done
            return seq_no_t

    return seq

def analyze_sequence(record, motifs_5prime, motifs_3prime):
    """
    Analyze a single sequence for required statistics.

    Args:
        record: BioPython SeqRecord object
        motifs_5prime: List of 5' motifs to check
        motifs_3prime: List of 3' motifs to check

    Returns:
        Dictionary with analysis results
    """
    sequence = str(record.seq).upper()
    seq_length = len(sequence)

    results = {
        'name': record.id,
        'has_5prime_motif': 0,
        'has_3prime_motif': 0,
        'longest_polyt': 0,
        'trailing_t_count': 0,
        'motif_5prime_type': 'none',  # Track which 5' motif was found
        'motif_3prime_type': 'none'   # Track which 3' motif was found
    }

    # 1. Check for 5' motifs within first 10 nt
    first_10nt = sequence[:min(10, seq_length)]
    motif_end = 0
    for motif in motifs_5prime:
        if motif in first_10nt:
            results['has_5prime_motif'] = 1
            results['motif_5prime_type'] = motif
            # Find the end position of the motif
            motif_start = first_10nt.index(motif)
            motif_end = motif_start + len(motif)
            break

    # 2. Count trailing Ts
    trailing_t_count = count_trailing_ts(sequence)
    results['trailing_t_count'] = trailing_t_count

    # 3. Properly trim polyT (including ...ATCGCATTTTTT cases) and check for 3' motif within last 12 nt of 3' end
    sequence_no_trailing_t = trim_polyt_proper(sequence)
    seq_no_t_length = len(sequence_no_trailing_t)

    # Add 1 T back to detect motifs ending in T (like ATCGT)
    seq_for_motif_check = sequence_no_trailing_t + 'T'

    # Check last 12 nt (or less if sequence is shorter)
    if seq_no_t_length > 0:
        last_12nt_start = max(0, len(seq_for_motif_check) - 12)
        last_12nt = seq_for_motif_check[last_12nt_start:]

        # Check for 3' motifs
        motif_start_in_seq = seq_no_t_length
        for motif in motifs_3prime:
            if motif in last_12nt:
                results['has_3prime_motif'] = 1
                results['motif_3prime_type'] = motif
                motif_pos_in_last12 = last_12nt.index(motif)
                motif_start_in_seq = last_12nt_start + motif_pos_in_last12
                break
    else:
        motif_start_in_seq = 0

    # 4. Find longest poly-T stretch between the two motifs
    # Region is from end of 5' motif to start of 3' motif
    if results['has_5prime_motif'] == 1 and results['has_3prime_motif'] == 1:
        # Both motifs present - search between them
        results['longest_polyt'] = find_longest_polyt_between_motifs(
            sequence, motif_end, motif_start_in_seq
        )
    elif results['has_5prime_motif'] == 1:
        # Only 5' motif present - search from end of 5' motif to end of sequence (minus trailing Ts)
        results['longest_polyt'] = find_longest_polyt_between_motifs(
            sequence, motif_end, seq_no_t_length
        )
    elif results['has_3prime_motif'] == 1:
        # Only 3' motif present - search from start to beginning of 3' motif
        results['longest_polyt'] = find_longest_polyt_between_motifs(
            sequence, 0, motif_start_in_seq
        )
    else:
        # No motifs present - search entire sequence (minus trailing Ts)
        results['longest_polyt'] = find_longest_polyt_between_motifs(
            sequence, 0, seq_no_t_length
        )

    return results

def process_fasta(input_source, output_handle, config_file='config.yaml'):
    """
    Process FASTA file and calculate statistics for each sequence.

    Args:
        input_source: Input file path or file handle
        output_handle: Output file handle
        config_file: Path to config file with motif definitions
    """
    # Load motifs from config
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        motifs_5prime = config.get('motifs', {}).get('5prime', ['GCGGT', 'GTGGT'])
        motifs_3prime = config.get('motifs', {}).get('3prime', ['ATCGC', 'ACCGC', 'ATCAC'])
    except Exception as e:
        print(f"Warning: Could not load config file, using defaults: {e}", file=sys.stderr)
        motifs_5prime = ['GCGGT', 'GTGGT']
        motifs_3prime = ['ATCGC', 'ACCGC', 'ATCAC']

    # Write header
    header = ['sequence_name', 'has_GCGGT_5prime', 'has_ATCGC_3prime',
              'longest_polyT_stretch', 'trailing_T_count', 'motif_5prime_type', 'motif_3prime_type']
    output_handle.write('\t'.join(header) + '\n')

    # Process sequences
    sequences_processed = 0

    # Handle input - could be file path or file handle
    if hasattr(input_source, 'read'):
        # It's a file handle
        for record in SeqIO.parse(input_source, 'fasta'):
            results = analyze_sequence(record, motifs_5prime, motifs_3prime)

            # Write results as TSV
            row = [
                results['name'],
                str(results['has_5prime_motif']),
                str(results['has_3prime_motif']),
                str(results['longest_polyt']),
                str(results['trailing_t_count']),
                results['motif_5prime_type'],
                results['motif_3prime_type']
            ]
            output_handle.write('\t'.join(row) + '\n')

            sequences_processed += 1
    else:
        # It's a file path
        with smart_open(input_source, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                results = analyze_sequence(record, motifs_5prime, motifs_3prime)

                # Write results as TSV
                row = [
                    results['name'],
                    str(results['has_5prime_motif']),
                    str(results['has_3prime_motif']),
                    str(results['longest_polyt']),
                    str(results['trailing_t_count']),
                    results['motif_5prime_type'],
                    results['motif_3prime_type']
                ]
                output_handle.write('\t'.join(row) + '\n')

                sequences_processed += 1

    print(f"Processed {sequences_processed} sequences", file=sys.stderr)

def main():
    parser = create_standard_parser(
        description='Calculate statistics for lncRNA sequences',
        epilog="""
Examples:
  # Process FASTA file
  python3 calculate_lncrna_stats.py lncrna.fa > lncrna.stats

  # Process from stdin
  cat lncrna.fa | python3 calculate_lncrna_stats.py > lncrna.stats

Output columns (TSV):
  1. sequence_name - Sequence identifier
  2. has_GCGGT_5prime - 1 if GCGGT or GTGGT found in first 10nt, 0 otherwise
  3. has_ATCGC_3prime - 1 if ATCGC found in last 10nt (after removing trailing Ts), 0 otherwise
  4. longest_polyT_stretch - Length of longest poly-T stretch between motifs
  5. trailing_T_count - Number of T nucleotides at 3' end
  6. motif_5prime_type - Type of 5' motif found: GCGGT, GTGGT, or none
        """
    )

    parser.add_argument('input', nargs='?',
                       help='Input FASTA file (use stdin if not provided)')
    parser.add_argument('-o', '--output',
                       help='Output TSV file (use stdout if not provided)')
    parser.add_argument('-c', '--config', default='config.yaml',
                       help='Configuration file with motif definitions (default: config.yaml)')

    args = parser.parse_args()

    # Determine input source
    if args.input:
        input_source = args.input
    else:
        input_source = sys.stdin

    # Determine output destination
    if args.output:
        output_handle = open(args.output, 'w')
    else:
        output_handle = sys.stdout

    try:
        process_fasta(input_source, output_handle, args.config)
    finally:
        if args.output:
            output_handle.close()

if __name__ == "__main__":
    main()