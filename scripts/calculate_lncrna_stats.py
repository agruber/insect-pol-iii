#!/usr/bin/env python3
"""
Calculate statistics for lncRNA sequences including motif presence and poly-T stretches.

Usage:
    python3 calculate_lncrna_stats.py input.fa > output.stats
    cat input.fa | python3 calculate_lncrna_stats.py > output.stats
"""

import sys
import re
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

def analyze_sequence(record):
    """
    Analyze a single sequence for required statistics.

    Args:
        record: BioPython SeqRecord object

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
        'trailing_t_count': 0
    }

    # 1. Check for GCGGT within first 10 nt of 5' end
    first_10nt = sequence[:min(10, seq_length)]
    if 'GCGGT' in first_10nt:
        results['has_5prime_motif'] = 1
        # Find the end position of the motif
        motif_start = first_10nt.index('GCGGT')
        motif_end = motif_start + 5
    else:
        motif_end = 0

    # 2. Count trailing Ts
    trailing_t_count = count_trailing_ts(sequence)
    results['trailing_t_count'] = trailing_t_count

    # 3. Remove trailing Ts and check for ATCGC within last 10 nt of 3' end
    if trailing_t_count > 0:
        sequence_no_trailing_t = sequence[:-trailing_t_count]
    else:
        sequence_no_trailing_t = sequence

    seq_no_t_length = len(sequence_no_trailing_t)

    # Check last 10 nt (or less if sequence is shorter)
    if seq_no_t_length > 0:
        last_10nt_start = max(0, seq_no_t_length - 10)
        last_10nt = sequence_no_trailing_t[last_10nt_start:]

        if 'ATCGC' in last_10nt:
            results['has_3prime_motif'] = 1
            # Find the start position of the motif in the original sequence
            motif_pos_in_last10 = last_10nt.index('ATCGC')
            motif_start_in_seq = last_10nt_start + motif_pos_in_last10
        else:
            motif_start_in_seq = seq_no_t_length
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

def process_fasta(input_source, output_handle):
    """
    Process FASTA file and calculate statistics for each sequence.

    Args:
        input_source: Input file path or file handle
        output_handle: Output file handle
    """
    # Write header
    header = ['sequence_name', 'has_GCGGT_5prime', 'has_ATCGC_3prime',
              'longest_polyT_stretch', 'trailing_T_count']
    output_handle.write('\t'.join(header) + '\n')

    # Process sequences
    sequences_processed = 0

    # Handle input - could be file path or file handle
    if hasattr(input_source, 'read'):
        # It's a file handle
        for record in SeqIO.parse(input_source, 'fasta'):
            results = analyze_sequence(record)

            # Write results as TSV
            row = [
                results['name'],
                str(results['has_5prime_motif']),
                str(results['has_3prime_motif']),
                str(results['longest_polyt']),
                str(results['trailing_t_count'])
            ]
            output_handle.write('\t'.join(row) + '\n')

            sequences_processed += 1
    else:
        # It's a file path
        with smart_open(input_source, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                results = analyze_sequence(record)

                # Write results as TSV
                row = [
                    results['name'],
                    str(results['has_5prime_motif']),
                    str(results['has_3prime_motif']),
                    str(results['longest_polyt']),
                    str(results['trailing_t_count'])
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
  2. has_GCGGT_5prime - 1 if GCGGT found in first 10nt, 0 otherwise
  3. has_ATCGC_3prime - 1 if ATCGC found in last 10nt (after removing trailing Ts), 0 otherwise
  4. longest_polyT_stretch - Length of longest poly-T stretch between motifs
  5. trailing_T_count - Number of T nucleotides at 3' end
        """
    )

    parser.add_argument('input', nargs='?',
                       help='Input FASTA file (use stdin if not provided)')
    parser.add_argument('-o', '--output',
                       help='Output TSV file (use stdout if not provided)')

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
        process_fasta(input_source, output_handle)
    finally:
        if args.output:
            output_handle.close()

if __name__ == "__main__":
    main()