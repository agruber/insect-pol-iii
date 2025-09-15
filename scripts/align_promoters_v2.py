#!/usr/bin/env python3
"""
Simplified promoter alignment script using k-mer anchors.
Reads FASTA from stdin, writes aligned FASTA to stdout.
"""

import sys
import os
import argparse
import math
import yaml
from typing import List, Dict, Set, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import numpy as np


def expand_iupac(kmer: str) -> List[str]:
    """Expand IUPAC ambiguity codes in a k-mer to all possible sequences."""
    iupac_codes = {
        'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
        'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
        'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
        'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
    }

    expanded = ['']
    for char in kmer.upper():
        if char in iupac_codes:
            new_expanded = []
            for seq in expanded:
                for base in iupac_codes[char]:
                    new_expanded.append(seq + base)
            expanded = new_expanded
        else:
            expanded = [seq + char for seq in expanded]

    return expanded


def generate_hamming_distance_1_variants(kmer: str) -> Set[str]:
    """
    Generate all Hamming distance 1 variants of a k-mer.

    Args:
        kmer: Input k-mer sequence

    Returns:
        Set of k-mers with Hamming distance 1 from input
    """
    variants = set()
    bases = ['A', 'C', 'G', 'T']
    kmer_upper = kmer.upper()

    for i in range(len(kmer_upper)):
        original_base = kmer_upper[i]
        for base in bases:
            if base != original_base:
                variant = kmer_upper[:i] + base + kmer_upper[i+1:]
                variants.add(variant)

    return variants


def find_kmer_in_region(sequence: str, kmer: str, search_start: int, search_end: int) -> Optional[int]:
    """
    Find a k-mer in the specified search region of a sequence.

    Args:
        sequence: DNA sequence to search
        kmer: k-mer to find
        search_start: Start position from end (e.g., -60 or -40)
        search_end: End position from end (e.g., -40 or -60)

    Returns:
        Position of k-mer from start of sequence if found, None otherwise
    """
    seq_len = len(sequence)
    kmer_len = len(kmer)

    # Handle both orderings of start/end
    # Ensure we search from the more negative to less negative position
    actual_start = min(search_start, search_end)  # More negative
    actual_end = max(search_start, search_end)    # Less negative

    # Convert negative positions to positive indices
    start_idx = max(0, seq_len + actual_start)
    end_idx = min(seq_len, seq_len + actual_end + kmer_len)  # Allow for k-mer length

    # Make sure the range is valid
    if start_idx >= end_idx:
        return None

    # Search for k-mer in the region
    search_region = sequence[start_idx:end_idx].upper()
    kmer_upper = kmer.upper()

    pos = search_region.find(kmer_upper)
    if pos != -1:
        return start_idx + pos

    return None


def find_anchor_positions(sequences: List[str], kmers: Set[str],
                         search_start: int, search_end: int) -> List[Optional[Tuple[int, str]]]:
    """
    Find anchor k-mer positions for each sequence.

    Returns:
        List of (position, kmer) tuples or None for each sequence
    """
    anchor_positions = []

    for seq in sequences:
        found = None
        for kmer in kmers:
            pos = find_kmer_in_region(seq, kmer, search_start, search_end)
            if pos is not None:
                found = (pos, kmer)
                break
        anchor_positions.append(found)

    return anchor_positions


def align_sequences_by_anchors(sequences: List[str],
                               anchor_data: List[Optional[Tuple[int, str]]]) -> List[str]:
    """
    Align sequences based on anchor positions.

    Args:
        sequences: Original sequences
        anchor_data: List of (position, kmer) tuples or None

    Returns:
        Aligned sequences with gaps
    """
    if not sequences:
        return []

    # Find the rightmost anchor position to determine alignment length
    max_right = 0
    max_left = 0

    for i, data in enumerate(anchor_data):
        if data is not None:
            pos, kmer = data
            seq_len = len(sequences[i])

            # Space needed on the left
            left_space = pos
            # Space needed on the right
            right_space = seq_len - pos

            max_left = max(max_left, left_space)
            max_right = max(max_right, right_space)

    # Total alignment length
    total_length = max_left + max_right

    # Create aligned sequences
    aligned = []
    for i, seq in enumerate(sequences):
        if anchor_data[i] is not None:
            pos, kmer = anchor_data[i]
            # Calculate padding
            left_pad = max_left - pos
            right_pad = total_length - left_pad - len(seq)

            # Create aligned sequence
            aligned_seq = '-' * left_pad + seq + '-' * right_pad
            aligned.append(aligned_seq)
        else:
            # No anchor - add as all gaps
            aligned.append('-' * total_length)

    return aligned


def mark_promoter_regions(aligned_seq: str, anchor_data: Optional[Tuple[int, str]],
                          promoter_shift: int, promoter_length: int, alignment_offset: int) -> str:
    """
    Mark promoter regions with uppercase, rest with lowercase.

    Args:
        aligned_seq: Aligned sequence with gaps
        anchor_data: (position, kmer) tuple or None
        promoter_shift: Nucleotides upstream of anchor where promoter starts
        promoter_length: Length of promoter region
        alignment_offset: Offset to convert original position to alignment position

    Returns:
        Sequence with promoter in uppercase, rest in lowercase
    """
    if anchor_data is None or '-' * len(aligned_seq) == aligned_seq:
        # No anchor or all gaps - return all lowercase
        return aligned_seq.lower()

    pos, kmer = anchor_data

    # Calculate promoter start position in alignment
    # pos is the position in original sequence
    # alignment_offset tells us how many gaps were added to the left
    anchor_pos_in_alignment = pos + alignment_offset
    promoter_start_in_alignment = anchor_pos_in_alignment - promoter_shift
    promoter_end_in_alignment = promoter_start_in_alignment + promoter_length

    # Ensure bounds are valid
    promoter_start_in_alignment = max(0, promoter_start_in_alignment)
    promoter_end_in_alignment = min(len(aligned_seq), promoter_end_in_alignment)

    # Build marked sequence
    marked = []
    for i, char in enumerate(aligned_seq):
        if char == '-':
            marked.append('-')
        elif promoter_start_in_alignment <= i < promoter_end_in_alignment:
            marked.append(char.upper())
        else:
            marked.append(char.lower())

    return ''.join(marked)


def load_rna_config(config_path: str = 'config.yaml') -> Dict:
    """Load RNA polymerase classification from config file."""
    # Try multiple possible paths for config.yaml
    possible_paths = [
        config_path,
        'config.yaml',
        os.path.join(os.path.dirname(__file__), '..', 'config.yaml'),
        '/home/andreas/WORK/SCIENCE/insect-pol-iii/config.yaml'
    ]

    for path in possible_paths:
        if os.path.exists(path):
            with open(path, 'r') as f:
                config = yaml.safe_load(f)
                return config

    # Return default config if file not found
    sys.stderr.write("Warning: config.yaml not found, using defaults\n")
    return {
        'rna_polymerases': {
            'pol2': {'types': ['U1', 'U2', 'U3', 'U4', 'U4atac', 'U5', 'U7', 'U11', 'U12']},
            'pol3': {'types': ['U6', 'U6atac', 'tRNAsec', 'RNase_MRP', 'Arthropod_7SK']}
        }
    }


def classify_sequence(seq_id: str, config: Dict) -> str:
    """Classify a sequence as Pol II or Pol III based on RNA type."""
    pol2_types = set(config['rna_polymerases']['pol2']['types'])
    pol3_types = set(config['rna_polymerases']['pol3']['types'])

    # Extract RNA type from sequence ID
    for rna_type in pol2_types:
        if f"|{rna_type}|" in seq_id:
            return 'pol2'

    for rna_type in pol3_types:
        if f"|{rna_type}|" in seq_id:
            return 'pol3'

    return 'unknown'


def build_pwm(sequences: List[str], promoter_start: int, promoter_length: int) -> Dict:
    """Build a Position Weight Matrix from aligned promoter sequences."""
    if not sequences:
        return None

    # Initialize frequency matrix
    freq_matrix = {'A': [], 'C': [], 'G': [], 'T': []}

    for pos in range(promoter_length):
        counts = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}  # Pseudocounts
        total = 0.4

        for seq in sequences:
            if promoter_start + pos < len(seq):
                nt = seq[promoter_start + pos].upper()
                if nt in 'ACGT':
                    counts[nt] += 1
                    total += 1

        # Normalize to frequencies
        for base in 'ACGT':
            freq_matrix[base].append(counts[base] / total if total > 0 else 0.25)

    # Calculate information content
    information = []
    for pos in range(promoter_length):
        ic = 0
        for base in 'ACGT':
            freq = freq_matrix[base][pos]
            if freq > 0:
                ic += freq * math.log2(freq / 0.25)
        information.append(max(0, ic))

    # Define core region (middle 60% of promoter)
    core_start = int(promoter_length * 0.2)
    core_end = int(promoter_length * 0.8)

    return {
        'length': promoter_length,
        'frequencies': freq_matrix,
        'information': information,
        'core_start': core_start,
        'core_end': core_end
    }


def calculate_match_score(sequence: str, pwm: Dict) -> float:
    """Calculate MATCH algorithm score for a sequence against a PWM."""
    if not pwm or not sequence:
        return 0.0

    seq_upper = sequence.upper()
    score = 0.0
    max_score = 0.0
    min_score = 0.0

    for i in range(min(len(seq_upper), pwm['length'])):
        if seq_upper[i] in 'ACGT':
            freq = pwm['frequencies'][seq_upper[i]][i]
            info = pwm['information'][i]

            # Information-weighted score
            if freq > 0:
                score += info * math.log2(4 * freq)
            else:
                score += -2 * info  # Penalty for zero frequency

            # Calculate max/min for this position
            max_freq = max(pwm['frequencies'][b][i] for b in 'ACGT')
            min_freq = min(pwm['frequencies'][b][i] for b in 'ACGT')

            if max_freq > 0:
                max_score += info * math.log2(4 * max_freq)
            if min_freq > 0:
                min_score += info * math.log2(4 * min_freq)
            else:
                min_score += -2 * info

    # Normalize to [0, 1]
    if max_score > min_score:
        normalized = (score - min_score) / (max_score - min_score)
        return max(0.0, min(1.0, normalized))

    return 0.5


def main():
    parser = argparse.ArgumentParser(
        description='Align promoter sequences using k-mer anchors',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--kmer', required=True,
                       help='K-mer(s) to use as anchors (supports IUPAC codes, comma-separated for multiple)')
    parser.add_argument('-s', '--search-start', type=int, default=-60,
                       help='Start of promoter search region from TSS')
    parser.add_argument('-e', '--search-end', type=int, default=-40,
                       help='End of promoter search region from TSS')
    parser.add_argument('-ps', '--promoter-shift', type=int, default=3,
                       help='Number of nucleotides upstream of anchor k-mer where promoter starts')
    parser.add_argument('-l', '--promoter-length', type=int, default=20,
                       help='Length of promoter region')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output to stderr')

    args = parser.parse_args()

    # Parse and expand k-mers
    initial_kmers = []
    expanded_kmers = set()

    for kmer in args.kmer.split(','):
        kmer = kmer.strip()
        initial_kmers.append(kmer)

        # Expand IUPAC codes
        expanded = expand_iupac(kmer)
        expanded_kmers.update(expanded)

        if args.verbose and len(expanded) > 1:
            sys.stderr.write(f"Expanded {kmer} -> {', '.join(expanded)}\n")

    # Always log the k-mers we're screening with
    sys.stderr.write(f"\nScreening with k-mers:\n")
    sys.stderr.write(f"  Initial: {', '.join(initial_kmers)}\n")
    if len(expanded_kmers) > 1:
        sys.stderr.write(f"  Expanded: {', '.join(sorted(expanded_kmers))}\n")
    sys.stderr.write(f"  Total: {len(expanded_kmers)} k-mer(s)\n")

    # Read sequences from stdin
    seq_records = list(SeqIO.parse(sys.stdin, 'fasta'))
    sequences = [str(record.seq) for record in seq_records]

    if args.verbose:
        sys.stderr.write(f"\nLoaded {len(sequences)} sequences\n")

    # Find anchor positions
    anchor_data = find_anchor_positions(sequences, expanded_kmers,
                                       args.search_start, args.search_end)

    # Count anchored sequences
    anchored_count = sum(1 for data in anchor_data if data is not None)

    sys.stderr.write(f"\nAnchoring results:\n")
    sys.stderr.write(f"  Anchored: {anchored_count}/{len(sequences)} sequences\n")

    # Always log k-mer positions for each sequence
    sys.stderr.write("\nK-mer positions:\n")
    for i, (record, data) in enumerate(zip(seq_records, anchor_data)):
        if data is not None:
            pos, kmer = data
            # Calculate position relative to end of sequence
            relative_pos = pos - len(sequences[i])
            sys.stderr.write(f"  {record.id}: {kmer}@{relative_pos}\n")
        else:
            sys.stderr.write(f"  {record.id}: none\n")

    # Generate HD1 variants for unanchored sequences
    unanchored_indices = [i for i, data in enumerate(anchor_data) if data is None]

    if unanchored_indices:
        # Generate all HD1 variants from the expanded k-mers
        hd1_kmers = set()
        for kmer in expanded_kmers:
            hd1_variants = generate_hamming_distance_1_variants(kmer)
            hd1_kmers.update(hd1_variants)

        # Remove original k-mers to get only new HD1 variants
        hd1_kmers -= expanded_kmers

        sys.stderr.write(f"\nGenerating HD1 variants for {len(unanchored_indices)} unanchored sequences:\n")
        sys.stderr.write(f"  Original k-mers: {len(expanded_kmers)}\n")
        sys.stderr.write(f"  HD1 variants: {len(hd1_kmers)}\n")
        sys.stderr.write(f"  Total k-mers for HD1 search: {len(expanded_kmers | hd1_kmers)}\n")

        # Find the best k-mer for each unanchored sequence
        best_hd1_anchors = {}  # {seq_index: (pos, kmer, score)}

        # We need PWMs to score, so do this after PWM building below

    # Align sequences
    aligned_sequences = align_sequences_by_anchors(sequences, anchor_data)

    # Load RNA classification config
    config = load_rna_config()

    # Classify sequences and separate by polymerase type
    pol2_sequences = []
    pol3_sequences = []
    pol2_aligned = []
    pol3_aligned = []

    for i, (record, data) in enumerate(zip(seq_records, anchor_data)):
        if data is not None:  # Only include anchored sequences
            seq_type = classify_sequence(record.id, config)
            if seq_type == 'pol2':
                pol2_sequences.append(record.id)
                pol2_aligned.append(aligned_sequences[i])
            elif seq_type == 'pol3':
                pol3_sequences.append(record.id)
                pol3_aligned.append(aligned_sequences[i])

    sys.stderr.write(f"\nSequence classification:\n")
    sys.stderr.write(f"  Pol II: {len(pol2_sequences)} sequences\n")
    sys.stderr.write(f"  Pol III: {len(pol3_sequences)} sequences\n")

    # Find promoter start position in alignment (same for all aligned sequences)
    promoter_start_in_alignment = None
    if anchor_data:
        for i, data in enumerate(anchor_data):
            if data is not None:
                pos, kmer = data
                # Count leading gaps
                left_pad = 0
                for char in aligned_sequences[i]:
                    if char == '-':
                        left_pad += 1
                    else:
                        break
                # Promoter starts at anchor position - promoter_shift
                promoter_start_in_alignment = left_pad + pos - args.promoter_shift
                break

    # Build PWMs from anchored sequences
    pol2_pwm = None
    pol3_pwm = None

    if pol2_aligned and promoter_start_in_alignment is not None:
        pol2_pwm = build_pwm(pol2_aligned, promoter_start_in_alignment, args.promoter_length)
        sys.stderr.write(f"  Built Pol II PWM from {len(pol2_aligned)} sequences\n")

    if pol3_aligned and promoter_start_in_alignment is not None:
        pol3_pwm = build_pwm(pol3_aligned, promoter_start_in_alignment, args.promoter_length)
        sys.stderr.write(f"  Built Pol III PWM from {len(pol3_aligned)} sequences\n")

    # Process unanchored sequences with HD1 k-mers
    if unanchored_indices and (pol2_pwm or pol3_pwm):
        all_hd1_kmers = expanded_kmers | hd1_kmers
        sys.stderr.write(f"\nProcessing HD1 k-mers for unanchored sequences:\n")

        for seq_idx in unanchored_indices:
            best_score = 0
            best_anchor = None
            record = seq_records[seq_idx]
            seq = sequences[seq_idx]

            # Test each HD1 k-mer
            for kmer in all_hd1_kmers:
                pos = find_kmer_in_region(seq, kmer, args.search_start, args.search_end)
                if pos is not None:
                    # Extract promoter region for this k-mer position
                    promoter_start = pos - args.promoter_shift
                    promoter_end = promoter_start + args.promoter_length

                    if promoter_start >= 0 and promoter_end <= len(seq):
                        promoter_region = seq[promoter_start:promoter_end]

                        # Score with both PWMs and take the better score
                        max_score = 0
                        if pol2_pwm:
                            pol2_score = calculate_match_score(promoter_region, pol2_pwm)
                            max_score = max(max_score, pol2_score)
                        if pol3_pwm:
                            pol3_score = calculate_match_score(promoter_region, pol3_pwm)
                            max_score = max(max_score, pol3_score)

                        if max_score > best_score:
                            best_score = max_score
                            best_anchor = (pos, kmer)

            # Only include if score > 0.8
            if best_score > 0.8:
                best_hd1_anchors[seq_idx] = (best_anchor[0], best_anchor[1], best_score)
                anchor_data[seq_idx] = best_anchor
                relative_pos = best_anchor[0] - len(seq)
                sys.stderr.write(f"  {record.id}: {best_anchor[1]}@{relative_pos} (score: {best_score:.3f})\n")
            else:
                sys.stderr.write(f"  {record.id}: no HD1 k-mer with score > 0.8 (best: {best_score:.3f})\n")

        # Re-align sequences with new anchors
        aligned_sequences = align_sequences_by_anchors(sequences, anchor_data)

        # Update anchored count
        new_anchored_count = sum(1 for data in anchor_data if data is not None)
        sys.stderr.write(f"  Updated anchored count: {new_anchored_count}/{len(sequences)} sequences\n")

        # Rebuild PWMs with the expanded set of anchored sequences
        pol2_sequences = []
        pol3_sequences = []
        pol2_aligned = []
        pol3_aligned = []

        for i, (record, data) in enumerate(zip(seq_records, anchor_data)):
            if data is not None:  # Only include anchored sequences
                seq_type = classify_sequence(record.id, config)
                if seq_type == 'pol2':
                    pol2_sequences.append(record.id)
                    pol2_aligned.append(aligned_sequences[i])
                elif seq_type == 'pol3':
                    pol3_sequences.append(record.id)
                    pol3_aligned.append(aligned_sequences[i])

        # Find new promoter start position in alignment
        promoter_start_in_alignment = None
        if anchor_data:
            for i, data in enumerate(anchor_data):
                if data is not None:
                    pos, kmer = data
                    # Count leading gaps
                    left_pad = 0
                    for char in aligned_sequences[i]:
                        if char == '-':
                            left_pad += 1
                        else:
                            break
                    # Promoter starts at anchor position - promoter_shift
                    promoter_start_in_alignment = left_pad + pos - args.promoter_shift
                    break

        # Rebuild PWMs with updated sequences
        if pol2_aligned and promoter_start_in_alignment is not None:
            pol2_pwm = build_pwm(pol2_aligned, promoter_start_in_alignment, args.promoter_length)
            sys.stderr.write(f"  Rebuilt Pol II PWM from {len(pol2_aligned)} sequences\n")

        if pol3_aligned and promoter_start_in_alignment is not None:
            pol3_pwm = build_pwm(pol3_aligned, promoter_start_in_alignment, args.promoter_length)
            sys.stderr.write(f"  Rebuilt Pol III PWM from {len(pol3_aligned)} sequences\n")

    # Calculate alignment offset (how many gaps were added to the left)
    # Find the maximum left padding
    max_left_pad = 0
    for i, data in enumerate(anchor_data):
        if data is not None:
            pos, kmer = data
            # Calculate left padding for this sequence
            left_pad = 0
            for char in aligned_sequences[i]:
                if char == '-':
                    left_pad += 1
                else:
                    break
            max_left_pad = max(max_left_pad, left_pad)

    # Score each sequence's promoter with both PWMs
    sys.stderr.write(f"\nPromoter scoring (MATCH algorithm):\n")
    sys.stderr.write(f"{'Sequence ID':<60} {'Type':<8} {'Pol2':<8} {'Pol3':<8}\n")
    sys.stderr.write("-" * 84 + "\n")

    # Determine alignment length for proper padding
    alignment_length = len(aligned_sequences[0]) if aligned_sequences else 0

    # Separate anchored and unanchored sequences
    anchored_outputs = []
    unanchored_outputs = []

    # Process all sequences and score anchored ones
    for i, (record, aligned_seq) in enumerate(zip(seq_records, aligned_sequences)):
        # Calculate alignment offset for this sequence
        left_pad = 0
        for char in aligned_seq:
            if char == '-':
                left_pad += 1
            else:
                break

        # Mark promoter regions
        marked_seq = mark_promoter_regions(
            aligned_seq, anchor_data[i],
            args.promoter_shift, args.promoter_length,
            left_pad
        )

        # Score promoter region if anchored
        if anchor_data[i] is not None and promoter_start_in_alignment is not None:
            # Extract promoter region for scoring
            promoter_end = promoter_start_in_alignment + args.promoter_length
            if promoter_start_in_alignment >= 0 and promoter_end <= len(aligned_seq):
                promoter_region = aligned_seq[promoter_start_in_alignment:promoter_end]

                # Calculate scores
                pol2_score = calculate_match_score(promoter_region, pol2_pwm) if pol2_pwm else 0.0
                pol3_score = calculate_match_score(promoter_region, pol3_pwm) if pol3_pwm else 0.0

                # Determine sequence type
                seq_type = classify_sequence(record.id, config)

                # Log scores
                sys.stderr.write(f"{record.id:<60} {seq_type:<8} {pol2_score:<8.3f} {pol3_score:<8.3f}\n")

            # Add to anchored outputs
            header = f">{record.id}"
            if record.description and record.description != record.id:
                header += f" {record.description}"
            anchored_outputs.append((header, marked_seq))

        else:
            # No anchor - no scoring
            sys.stderr.write(f"{record.id:<60} {'none':<8} {'--':<8} {'--':<8}\n")

            # Create unanchored sequence output (all lowercase, padded at 3' end)
            original_seq = str(seq_records[i].seq).lower()
            # Pad with gaps at 3' end to match alignment length
            gaps_needed = alignment_length - len(original_seq)
            padded_seq = original_seq + '-' * gaps_needed

            header = f">{record.id}"
            if record.description and record.description != record.id:
                header += f" {record.description}"
            unanchored_outputs.append((header, padded_seq))

    # Output anchored sequences first
    for header, sequence in anchored_outputs:
        sys.stdout.write(f"{header}\n{sequence}\n")

    # Output unanchored sequences at the end
    for header, sequence in unanchored_outputs:
        sys.stdout.write(f"{header}\n{sequence}\n")

    if args.verbose:
        sys.stderr.write(f"\nAlignment complete!\n")
        sys.stderr.write(f"  Total sequences: {len(aligned_sequences)}\n")
        sys.stderr.write(f"  Alignment length: {len(aligned_sequences[0]) if aligned_sequences else 0} bp\n")


if __name__ == '__main__':
    main()