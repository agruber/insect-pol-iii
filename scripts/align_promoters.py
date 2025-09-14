#!/usr/bin/env python3
"""
Align upstream sequences by anchoring conserved k-mers and detect promoter regions.

This script performs multiple sequence alignment (MSA) based on conserved k-mer anchors
in promoter regions, then identifies and annotates the most conserved promoter window.
"""

import sys
import argparse
from collections import defaultdict, Counter
from typing import List, Tuple, Dict, Set, Optional
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools


# IUPAC nucleotide codes
IUPAC_CODES = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'U': ['T'],  # Treat U as T
    'R': ['A', 'G'],  # puRine
    'Y': ['C', 'T'],  # pYrimidine
    'S': ['G', 'C'],  # Strong interaction
    'W': ['A', 'T'],  # Weak interaction
    'K': ['G', 'T'],  # Keto
    'M': ['A', 'C'],  # aMino
    'B': ['C', 'G', 'T'],  # not A
    'D': ['A', 'G', 'T'],  # not C
    'H': ['A', 'C', 'T'],  # not G
    'V': ['A', 'C', 'G'],  # not T
    'N': ['A', 'C', 'G', 'T']  # aNy
}


def expand_iupac_kmer(kmer: str) -> Set[str]:
    """
    Expand an IUPAC-encoded k-mer into all possible DNA sequences.

    Args:
        kmer: K-mer with IUPAC codes (e.g., 'RTTT' -> R can be A or G)

    Returns:
        Set of all possible DNA sequences

    Examples:
        'RTTT' -> {'ATTT', 'GTTT'}
        'TATAWW' -> {'TATAAA', 'TATAAT', 'TATATA', 'TATATT'}
    """
    kmer = kmer.upper()

    # Get all possible nucleotides for each position
    position_options = []
    for char in kmer:
        if char in IUPAC_CODES:
            position_options.append(IUPAC_CODES[char])
        else:
            # Unknown character, treat as N (any nucleotide)
            sys.stderr.write(f"Warning: Unknown IUPAC code '{char}', treating as N\n")
            position_options.append(['A', 'C', 'G', 'T'])

    # Generate all combinations
    expanded = set()
    for combo in itertools.product(*position_options):
        expanded.add(''.join(combo))

    return expanded


def extract_kmers_from_region(sequence: str, k: int, start: int, end: int) -> Dict[str, List[int]]:
    """
    Extract k-mers from a specific region of the sequence.
    
    Args:
        sequence: DNA sequence
        k: k-mer length
        start: Start position (negative, from end)
        end: End position (negative, from end)
    
    Returns:
        Dictionary mapping k-mer to list of positions where it occurs
    """
    kmers = defaultdict(list)
    seq_len = len(sequence)
    
    # Convert negative positions to positive indices
    region_start = max(0, seq_len + start)
    region_end = min(seq_len, seq_len + end)
    
    if region_start >= region_end or region_end - region_start < k:
        return kmers
    
    for i in range(region_start, region_end - k + 1):
        kmer = sequence[i:i+k].upper()
        # Store position relative to sequence start
        kmers[kmer].append(i)
    
    return kmers


def score_kmers(sequences: List[str], k: int, start: int, end: int) -> List[Tuple[str, float, int, float]]:
    """
    Score k-mers by coverage and positional conservation.
    
    Args:
        sequences: List of DNA sequences
        k: k-mer length
        start: Start of search region (negative)
        end: End of search region (negative)
    
    Returns:
        List of (kmer, score, coverage, position_std) sorted by score
    """
    kmer_positions = defaultdict(list)
    kmer_sequences = defaultdict(set)
    
    for seq_idx, seq in enumerate(sequences):
        kmers = extract_kmers_from_region(seq, k, start, end)
        for kmer, positions in kmers.items():
            # Track which sequences contain this k-mer
            kmer_sequences[kmer].add(seq_idx)
            # Track all positions (relative to end of sequence)
            for pos in positions:
                kmer_positions[kmer].append(pos - len(seq))
    
    scored_kmers = []
    for kmer in kmer_positions:
        coverage = len(kmer_sequences[kmer])
        positions = kmer_positions[kmer]
        
        # Calculate positional conservation (lower std = more conserved)
        if len(positions) > 1:
            position_std = np.std(positions)
        else:
            position_std = 0.0
        
        # Score: high coverage, low positional variance
        # Normalize std by region size to make it comparable
        region_size = abs(end - start)
        normalized_std = position_std / region_size if region_size > 0 else 0
        score = coverage / len(sequences) - normalized_std * 0.5
        
        scored_kmers.append((kmer, score, coverage, position_std))
    
    return sorted(scored_kmers, key=lambda x: x[1], reverse=True)


def hamming_distance(s1: str, s2: str) -> int:
    """
    Calculate Hamming distance between two strings of equal length.
    
    Args:
        s1: First string
        s2: Second string
    
    Returns:
        Number of positions where characters differ
    """
    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def find_related_kmers(kmer: str, max_mismatches: int = 1) -> Set[str]:
    """
    Find k-mers with up to max_mismatches differences.
    
    Args:
        kmer: Original k-mer
        max_mismatches: Maximum number of mismatches allowed
    
    Returns:
        Set of related k-mers
    """
    related = {kmer}
    nucleotides = 'ACGT'
    
    if max_mismatches >= 1:
        # Generate 1-mismatch variants
        for i in range(len(kmer)):
            for nt in nucleotides:
                if nt != kmer[i]:
                    variant = kmer[:i] + nt + kmer[i+1:]
                    related.add(variant)
    
    return related


def find_anchor_positions(sequences: List[str], anchor_kmers: Set[str], k: int) -> List[Optional[int]]:
    """
    Find positions of anchor k-mers in each sequence.
    
    Args:
        sequences: List of DNA sequences
        anchor_kmers: Set of k-mers to search for
        k: k-mer length
    
    Returns:
        List of positions (or None if not found) for each sequence
    """
    positions = []
    
    for seq in sequences:
        seq_upper = seq.upper()
        best_pos = None
        
        # Search for any anchor k-mer
        for kmer in anchor_kmers:
            pos = seq_upper.find(kmer)
            if pos != -1:
                if best_pos is None or pos < best_pos:
                    best_pos = pos
        
        positions.append(best_pos)
    
    return positions


def align_sequences_by_anchors(sequences: List[str], anchor_positions: List[Optional[int]], k: int) -> Tuple[List[str], List[int], List[int]]:
    """
    Align sequences by placing anchor k-mers at the same column.
    Sequences without anchors are kept separate.
    
    Args:
        sequences: List of DNA sequences
        anchor_positions: Position of anchor in each sequence
        k: k-mer length
    
    Returns:
        Tuple of:
        - List of all aligned sequences (anchored first, then non-anchored)
        - List of indices for anchored sequences
        - List of indices for non-anchored sequences
    """
    # Separate sequences with and without anchors
    anchored_indices = [i for i, pos in enumerate(anchor_positions) if pos is not None]
    non_anchored_indices = [i for i, pos in enumerate(anchor_positions) if pos is None]
    
    if not anchored_indices:
        # No anchors found, return original sequences with gaps
        max_len = max(len(s) for s in sequences) if sequences else 0
        aligned = [s + '-' * (max_len - len(s)) for s in sequences]
        return aligned, [], list(range(len(sequences)))
    
    # Find the maximum anchor position to use as alignment target
    max_anchor_pos = max(anchor_positions[i] for i in anchored_indices)
    
    # First, align sequences WITH anchors
    aligned_anchored = []
    for i in anchored_indices:
        # Add gaps to align anchor at max_anchor_pos position
        left_pad = max_anchor_pos - anchor_positions[i]
        aligned_seq = '-' * left_pad + sequences[i]
        aligned_anchored.append(aligned_seq)
    
    # Find max length of anchored sequences
    max_len = max(len(s) for s in aligned_anchored) if aligned_anchored else 0
    
    # Pad anchored sequences to same length
    aligned_anchored = [s + '-' * (max_len - len(s)) for s in aligned_anchored]
    
    # Now handle sequences WITHOUT anchors - include sequence and pad with gaps
    aligned_non_anchored = []
    for i in non_anchored_indices:
        # Add the original sequence and pad with gaps to match alignment length
        seq_len = len(sequences[i])
        if seq_len < max_len:
            # Pad sequence to match alignment length
            aligned_seq = sequences[i] + '-' * (max_len - seq_len)
        else:
            # If sequence is longer than current max, use it as is
            # and we'll need to update max_len
            aligned_seq = sequences[i]
        aligned_non_anchored.append(aligned_seq)
    
    # Update max_len if any non-anchored sequence is longer
    if aligned_non_anchored:
        max_len = max(max_len, max(len(s) for s in aligned_non_anchored))
        # Re-pad all sequences to new max length
        aligned_anchored = [s + '-' * (max_len - len(s)) for s in aligned_anchored]
        aligned_non_anchored = [s + '-' * (max_len - len(s)) for s in aligned_non_anchored]
    
    # Combine all sequences: anchored first, then non-anchored
    all_aligned = aligned_anchored + aligned_non_anchored
    
    return all_aligned, anchored_indices, non_anchored_indices


def calculate_conservation_score(aligned_seqs: List[str], window_start: int, window_end: int) -> float:
    """
    Calculate conservation score for a window in the alignment.
    
    Args:
        aligned_seqs: List of aligned sequences
        window_start: Start position of window
        window_end: End position of window
    
    Returns:
        Conservation score (0-1)
    """
    if not aligned_seqs or window_start >= window_end:
        return 0.0
    
    scores = []
    for pos in range(window_start, min(window_end, len(aligned_seqs[0]))):
        # Count nucleotides at this position
        counts = Counter()
        total = 0
        for seq in aligned_seqs:
            if pos < len(seq) and seq[pos] != '-':
                counts[seq[pos].upper()] += 1
                total += 1
        
        if total > 0:
            # Shannon entropy-based conservation
            max_count = max(counts.values())
            conservation = max_count / total
            scores.append(conservation)
    
    return np.mean(scores) if scores else 0.0


def find_best_promoter_window(aligned_seqs: List[str], anchor_positions: List[Optional[int]], 
                              window_size: int, k: int, search_start: int, search_end: int,
                              original_sequences: List[str]) -> Tuple[int, int]:
    """
    Find the best conserved promoter window within the search region.
    
    Args:
        aligned_seqs: List of aligned sequences
        anchor_positions: Original anchor positions in unaligned sequences
        window_size: Size of promoter window
        k: k-mer length
        search_start: Start of search region (negative, from end)
        search_end: End of search region (negative, from end)
        original_sequences: Original unaligned sequences
    
    Returns:
        Tuple of (start, end) positions for best window in aligned sequences
    """
    if not aligned_seqs:
        return (0, 0)
    
    # Calculate the search region boundaries in the aligned sequences
    # We need to map the original search region to aligned coordinates
    valid_indices = [i for i, pos in enumerate(anchor_positions) if pos is not None]
    if not valid_indices:
        # No anchors found, can't properly constrain the search
        return (0, min(window_size, len(aligned_seqs[0])))
    
    # For each aligned sequence, find where the search region maps to
    search_regions_aligned = []
    for i, aligned_seq in enumerate(aligned_seqs):
        if i in valid_indices:
            # This sequence has an anchor
            orig_len = len(original_sequences[i])
            # Calculate original search region positions
            orig_search_start = max(0, orig_len + search_start)
            orig_search_end = min(orig_len, orig_len + search_end)
            
            # Find these positions in the aligned sequence
            # Count non-gap positions to map back
            orig_pos = 0
            aligned_start = None
            aligned_end = None
            
            for j, char in enumerate(aligned_seq):
                if char != '-':
                    if orig_pos == orig_search_start:
                        aligned_start = j
                    if orig_pos == orig_search_end:
                        aligned_end = j
                        break
                    orig_pos += 1
            
            if aligned_start is not None and aligned_end is not None:
                search_regions_aligned.append((aligned_start, aligned_end))
    
    if not search_regions_aligned:
        # Fallback to anchor-based search
        anchor_col = max(anchor_positions[i] for i in valid_indices)
        best_window = (max(0, anchor_col - window_size // 2), 
                       min(len(aligned_seqs[0]), anchor_col + window_size // 2))
        return best_window
    
    # Find the common search region across all aligned sequences
    min_start = min(s for s, e in search_regions_aligned)
    max_end = max(e for s, e in search_regions_aligned)
    
    # Search for best window within these constraints
    best_score = 0
    best_window = (min_start, min(min_start + window_size, max_end))
    
    # Ensure window size doesn't exceed the search region
    actual_window_size = min(window_size, max_end - min_start)
    
    for start in range(min_start, max(min_start + 1, max_end - actual_window_size + 1)):
        end = min(start + actual_window_size, max_end)
        score = calculate_conservation_score(aligned_seqs, start, end)
        if score > best_score:
            best_score = score
            best_window = (start, end)
    
    return best_window


def annotate_promoter_regions(aligned_seqs: List[str], promoter_start: int, promoter_end: int) -> List[str]:
    """
    Annotate promoter regions with uppercase, rest with lowercase.
    
    Args:
        aligned_seqs: List of aligned sequences
        promoter_start: Start of promoter region
        promoter_end: End of promoter region
    
    Returns:
        List of annotated sequences
    """
    annotated = []
    
    for seq in aligned_seqs:
        if promoter_start >= len(seq) or promoter_end <= 0:
            annotated.append(seq.lower())
            continue
        
        # Split sequence into regions
        before = seq[:promoter_start].lower() if promoter_start > 0 else ''
        promoter = seq[promoter_start:promoter_end].upper()
        after = seq[promoter_end:].lower() if promoter_end < len(seq) else ''
        
        annotated.append(before + promoter + after)
    
    return annotated


def main():
    parser = argparse.ArgumentParser(
        description='Align upstream sequences by conserved k-mer anchors and detect promoters',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Automatic k-mer discovery
  cat upstream.fasta | python3 align_promoters.py > aligned.fasta
  
  # Custom k-mer anchors with IUPAC codes
  cat upstream.fasta | python3 align_promoters.py --kmer TATAAA --kmer TATAWW > aligned.fasta

  # Example: RTTT expands to ATTT and GTTT (R = A or G)
  cat upstream.fasta | python3 align_promoters.py --kmer RTTT > aligned.fasta
  
  # Custom parameters
  cat upstream.fasta | python3 align_promoters.py -k 6 -s -45 -e -55 -l 25 > aligned.fasta
        '''
    )
    
    parser.add_argument('-k', '--kmer-size', type=int, default=6,
                        help='K-mer size for anchor discovery (default: 6)')
    parser.add_argument('-s', '--search-start', type=int, default=-40,
                        help='Start of promoter search region from TSS (default: -40)')
    parser.add_argument('-e', '--search-end', type=int, default=-60,
                        help='End of promoter search region from TSS (default: -60)')
    parser.add_argument('-l', '--promoter-length', type=int, default=20,
                        help='Length of promoter window (default: 20)')
    parser.add_argument('--kmer', action='append', dest='kmers',
                        help='Specific k-mer to use as anchor (supports IUPAC codes, can be specified multiple times)')
    parser.add_argument('-m', '--min-coverage', type=int, default=3,
                        help='Minimum sequences that must contain k-mer (default: 3)')
    parser.add_argument('-t', '--top', type=int, default=1,
                        help='Use top N k-mers that have 1 nt distance to the best k-mer (default: 1)')
    parser.add_argument('-o', '--output', type=str, default='-',
                        help='Output file (default: stdout)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress to stderr')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.search_start > args.search_end:
        args.search_start, args.search_end = args.search_end, args.search_start
    
    if args.promoter_length < args.kmer_size:
        sys.stderr.write(f"Error: Promoter length ({args.promoter_length}) must be >= k-mer size ({args.kmer_size})\n")
        sys.exit(1)
    
    # Read sequences from stdin
    sequences = []
    seq_records = []
    for record in SeqIO.parse(sys.stdin, 'fasta'):
        sequences.append(str(record.seq))
        seq_records.append(record)
    
    if not sequences:
        sys.stderr.write("Error: No sequences provided\n")
        sys.exit(1)
    
    if args.verbose:
        sys.stderr.write(f"Loaded {len(sequences)} sequences\n")
    
    # Step 1-3: K-mer discovery or use provided k-mers
    if args.kmers:
        # Expand IUPAC codes in provided k-mers
        anchor_kmers = set()
        sys.stderr.write(f"\nProvided k-mers with IUPAC codes:\n")
        for kmer in args.kmers:
            expanded = expand_iupac_kmer(kmer)
            anchor_kmers.update(expanded)
            sys.stderr.write(f"  {kmer} -> {', '.join(sorted(expanded))}\n")
        sys.stderr.write(f"\nTotal expanded k-mers: {len(anchor_kmers)}\n")
        sys.stderr.write(f"Using k-mers: {', '.join(sorted(anchor_kmers))}\n")
    else:
        # Discover k-mers
        sys.stderr.write(f"\nSearching for {args.kmer_size}-mers in region [{args.search_start}, {args.search_end}] relative to sequence end...\n")
        scored_kmers = score_kmers(sequences, args.kmer_size, args.search_start, args.search_end)
        
        if not scored_kmers:
            sys.stderr.write("Error: No k-mers found in search region\n")
            sys.exit(1)
        
        # Show top k-mers
        sys.stderr.write(f"\nTop {min(10, len(scored_kmers))} k-mers by score:\n")
        sys.stderr.write("K-mer\tScore\tCoverage\tPos.Std\n")
        for i, (kmer, score, coverage, pos_std) in enumerate(scored_kmers[:10]):
            sys.stderr.write(f"{kmer}\t{score:.3f}\t{coverage}/{len(sequences)}\t{pos_std:.2f}\n")
        
        # Select best k-mer that meets coverage threshold
        best_kmer = None
        for kmer, score, coverage, pos_std in scored_kmers:
            if coverage >= args.min_coverage:
                best_kmer = kmer
                sys.stderr.write(f"\nSelected best k-mer: {kmer}\n")
                sys.stderr.write(f"  Score: {score:.3f}\n")
                sys.stderr.write(f"  Coverage: {coverage}/{len(sequences)} sequences\n")
                sys.stderr.write(f"  Position StdDev: {pos_std:.2f} bp\n")
                break
        
        if not best_kmer:
            sys.stderr.write(f"\nError: No k-mer found with minimum coverage of {args.min_coverage} sequences\n")
            sys.stderr.write(f"Best k-mer had coverage of {scored_kmers[0][2]} sequences\n")
            sys.exit(1)
        
        # Collect top N k-mers that have 1 nt distance to the best k-mer
        anchor_kmers = {best_kmer}
        selected_kmers = [best_kmer]
        
        if args.top > 1:
            sys.stderr.write(f"\nSelecting top {args.top} k-mers with Hamming distance ≤1 to best k-mer:\n")
            candidates_added = 0
            for kmer, score, coverage, pos_std in scored_kmers:
                if kmer == best_kmer:
                    continue
                if hamming_distance(kmer, best_kmer) == 1:
                    if coverage >= args.min_coverage:
                        anchor_kmers.add(kmer)
                        selected_kmers.append(kmer)
                        candidates_added += 1
                        sys.stderr.write(f"  Added: {kmer} (score: {score:.3f}, coverage: {coverage}/{len(sequences)})\n")
                        if len(selected_kmers) >= args.top:
                            break
            
            if candidates_added == 0:
                sys.stderr.write("  No additional k-mers with Hamming distance 1 found\n")
        
        # Use only the selected k-mers as anchors (no expansion)
        anchor_kmers = set(selected_kmers)
        sys.stderr.write(f"\nUsing {len(selected_kmers)} top k-mer(s): {', '.join(selected_kmers)}\n")
    
    # Step 4: Find anchor positions and align
    anchor_positions = find_anchor_positions(sequences, anchor_kmers, args.kmer_size)
    
    # Count how many sequences have anchors and show positions
    anchored_count = sum(1 for pos in anchor_positions if pos is not None)
    sys.stderr.write(f"\nAnchor k-mers found in {anchored_count}/{len(sequences)} sequences\n")
    
    if anchored_count == 0:
        sys.stderr.write("Error: No anchor k-mers found in sequences\n")
        sys.exit(1)
    
    if args.verbose:
        sys.stderr.write("\nAnchor positions in each sequence:\n")
        for i, (record, pos) in enumerate(zip(seq_records, anchor_positions)):
            if pos is not None:
                relative_pos = pos - len(sequences[i])
                sys.stderr.write(f"  {record.id}: position {pos} (relative: {relative_pos})\n")
            else:
                sys.stderr.write(f"  {record.id}: no anchor found\n")
    
    sys.stderr.write(f"\nAligning sequences by anchor positions...\n")
    aligned_sequences, anchored_indices, non_anchored_indices = align_sequences_by_anchors(
        sequences, anchor_positions, args.kmer_size
    )
    
    # Report if any sequences lack anchors
    if non_anchored_indices:
        sys.stderr.write(f"  {len(non_anchored_indices)} sequence(s) without anchors will be appended with gaps only\n")
        if args.verbose:
            for idx in non_anchored_indices:
                sys.stderr.write(f"    - {seq_records[idx].id}\n")
    
    # Step 5: Find best promoter window (only for anchored sequences)
    if anchored_indices:
        sys.stderr.write(f"Searching for best conserved {args.promoter_length}bp promoter window within search region [{args.search_start}, {args.search_end}]...\n")
        
        # Only use anchored sequences for finding promoter window
        anchored_seqs_aligned = aligned_sequences[:len(anchored_indices)]
        anchored_positions = [anchor_positions[i] for i in anchored_indices]
        anchored_seqs_original = [sequences[i] for i in anchored_indices]
        
        promoter_start, promoter_end = find_best_promoter_window(
            anchored_seqs_aligned, anchored_positions, args.promoter_length, args.kmer_size,
            args.search_start, args.search_end, anchored_seqs_original
        )
        
        conservation = calculate_conservation_score(anchored_seqs_aligned, promoter_start, promoter_end)
        sys.stderr.write(f"\nBest promoter window: columns {promoter_start}-{promoter_end}\n")
        sys.stderr.write(f"  Conservation score: {conservation:.3f}\n")
        sys.stderr.write(f"  Window size: {promoter_end - promoter_start} bp\n")
    else:
        # No anchored sequences, no promoter to annotate
        promoter_start, promoter_end = 0, 0
        sys.stderr.write("\nNo anchored sequences found - no promoter window to annotate\n")
    
    # Step 6: Annotate promoter regions (only for anchored sequences)
    annotated_sequences = []
    
    # Annotate anchored sequences with promoter regions
    for i in range(len(anchored_indices)):
        if promoter_start < promoter_end:
            annotated = annotate_promoter_regions([aligned_sequences[i]], promoter_start, promoter_end)[0]
        else:
            annotated = aligned_sequences[i].lower()
        annotated_sequences.append(annotated)
    
    # Add non-anchored sequences without annotation (just lowercase/gaps)
    for i in range(len(anchored_indices), len(aligned_sequences)):
        annotated_sequences.append(aligned_sequences[i].lower())
    
    # Output aligned and annotated sequences
    out_handle = sys.stdout if args.output == '-' else open(args.output, 'w')
    
    # Output anchored sequences first
    for i, idx in enumerate(anchored_indices):
        record = seq_records[idx]
        annotated_seq = annotated_sequences[i]
        # Write FASTA format directly to avoid line wrapping
        out_handle.write(f">{record.id}")
        if record.description and record.description != record.id:
            out_handle.write(f" {record.description}")
        out_handle.write(f"\n{annotated_seq}\n")
    
    # Output non-anchored sequences
    for i, idx in enumerate(non_anchored_indices):
        record = seq_records[idx]
        annotated_seq = annotated_sequences[len(anchored_indices) + i]
        # Write FASTA format directly to avoid line wrapping
        out_handle.write(f">{record.id}")
        if record.description and record.description != record.id:
            out_handle.write(f" {record.description}")
        out_handle.write(f" [no anchor]\n{annotated_seq}\n")
    
    if args.output != '-':
        out_handle.close()
    
    sys.stderr.write(f"\nAlignment complete!\n")
    sys.stderr.write(f"  Total sequences: {len(annotated_sequences)}\n")
    sys.stderr.write(f"    - Anchored (with promoter): {len(anchored_indices)}\n")
    sys.stderr.write(f"    - Non-anchored (gaps only): {len(non_anchored_indices)}\n")
    sys.stderr.write(f"  Alignment length: {len(annotated_sequences[0]) if annotated_sequences else 0} bp\n")
    sys.stderr.write(f"  Output written to: {args.output if args.output != '-' else 'stdout'}\n")


if __name__ == '__main__':
    main()