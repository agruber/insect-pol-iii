#!/usr/bin/env python3
"""
K-mer based promoter region detection and MSA generation.
Uses position-specific k-mer analysis to find conserved promoter elements
and builds MSA with promoter regions in uppercase.

Algorithm:
1. Extract k-mers from promoter region (-60 to -40 relative to TSS)
2. Find most conserved k-mer as anchor point
3. Align sequences using anchor-based positioning
4. Generate MSA with promoter regions (length=20) in uppercase

Usage:
    python3 annotate_promoter_regions.py -i upstream.fa -o promoters_msa.fa
"""

import sys
import argparse
from Bio import SeqIO
from collections import defaultdict, Counter
import itertools

def extract_promoter_kmers(sequences, k=6, promoter_start=-60, promoter_end=-40):
    """Extract k-mers from promoter region relative to TSS (sequence end)"""
    kmer_positions = defaultdict(list)  # kmer -> list of (seq_idx, position_in_promoter_region)
    
    for seq_idx, seq in enumerate(sequences):
        seq_str = str(seq).upper()
        seq_len = len(seq_str)
        
        # Convert relative positions to absolute positions
        # TSS is at sequence end, so -60 to -40 means positions (seq_len-60) to (seq_len-40)
        abs_start = seq_len + promoter_start  # promoter_start is negative
        abs_end = seq_len + promoter_end
        
        # Skip sequences too short for promoter region
        if abs_start < 0 or abs_end > seq_len or abs_start >= abs_end:
            continue
        
        promoter_seq = seq_str[abs_start:abs_end]
        
        # Extract all k-mers from this promoter region
        for i in range(len(promoter_seq) - k + 1):
            kmer = promoter_seq[i:i+k]
            if 'N' not in kmer:  # Skip ambiguous k-mers
                relative_pos = i  # Position within promoter region
                kmer_positions[kmer].append((seq_idx, relative_pos))
    
    return kmer_positions

def hamming_distance(seq1, seq2):
    """Calculate Hamming distance between two sequences of equal length"""
    if len(seq1) != len(seq2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def find_best_anchor_kmers(kmer_positions, min_sequences=3, top_kmers=3):
    """Find the best anchor k-mer and related k-mers with 1 nucleotide difference"""
    kmer_scores = []
    
    # First, find all qualified k-mers and their scores
    for kmer, positions in kmer_positions.items():
        # Count how many sequences contain this k-mer
        seq_count = len(set(pos[0] for pos in positions))
        
        if seq_count >= min_sequences:
            # Calculate position consistency (how clustered are the positions)
            relative_positions = [pos[1] for pos in positions]
            pos_counter = Counter(relative_positions)
            
            # Score based on sequence coverage and position consistency
            max_pos_count = max(pos_counter.values())
            consistency_score = max_pos_count / len(positions)
            coverage_score = seq_count
            
            total_score = coverage_score * consistency_score
            
            kmer_scores.append((kmer, total_score, seq_count))
    
    if not kmer_scores:
        return []
    
    # Sort by score (best first)
    kmer_scores.sort(key=lambda x: x[1], reverse=True)
    
    # Get the best k-mer as the primary anchor
    best_kmer = kmer_scores[0][0]
    selected_kmers = [kmer_scores[0]]  # Start with the best k-mer
    
    # Find additional k-mers that have exactly 1 nucleotide difference from the best
    for kmer, score, coverage in kmer_scores[1:]:
        if len(selected_kmers) >= top_kmers:
            break
            
        # Check if this k-mer has exactly 1 nucleotide difference from the best
        if hamming_distance(best_kmer, kmer) == 1:
            selected_kmers.append((kmer, score, coverage))
    
    return selected_kmers

def align_sequences_with_multiple_anchors(sequences, headers, top_anchor_kmers, kmer_positions, 
                                        promoter_length=20, max_gaps=40):
    """Align sequences using multiple anchor k-mers, ensuring all anchors are at the same MSA position"""
    
    # First, find the global target position for ALL anchors
    all_anchor_positions = []
    sequence_anchor_map = {}  # seq_idx -> (anchor_kmer, anchor_pos_absolute)
    
    # Collect all anchor positions from all sequences
    for anchor_kmer, _, _ in top_anchor_kmers:
        for seq_idx, relative_pos in kmer_positions[anchor_kmer]:
            seq = sequences[seq_idx]
            seq_str = str(seq).upper()
            seq_len = len(seq_str)
            
            # Find absolute position of k-mer in the full sequence
            promoter_start_abs = seq_len - 60
            anchor_pos_absolute = promoter_start_abs + relative_pos
            
            all_anchor_positions.append(anchor_pos_absolute)
            
            # Store best anchor for each sequence (first match wins = highest priority)
            if seq_idx not in sequence_anchor_map:
                sequence_anchor_map[seq_idx] = (anchor_kmer, anchor_pos_absolute)
    
    if not all_anchor_positions:
        return [], [], [], set()
    
    # Use MAXIMUM position as global target so all anchors align at same MSA position
    global_target_pos = max(all_anchor_positions)
    
    print(f"Global anchor target position: {global_target_pos}", file=sys.stderr)
    for anchor_kmer, _, _ in top_anchor_kmers:
        count = sum(1 for seq_idx, _ in kmer_positions[anchor_kmer])
        print(f"  {anchor_kmer}: {count} sequences", file=sys.stderr)
    
    # Process all sequences with anchors
    aligned_sequences = []
    aligned_headers = []
    alignment_scores = []
    used_sequences = set()
    
    # Group sequences by their assigned anchor
    anchor_groups = {}
    for seq_idx, (anchor_kmer, anchor_pos_absolute) in sequence_anchor_map.items():
        if anchor_kmer not in anchor_groups:
            anchor_groups[anchor_kmer] = []
        anchor_groups[anchor_kmer].append((seq_idx, anchor_pos_absolute))
    
    # Process each anchor group
    for anchor_rank, (anchor_kmer, anchor_score, anchor_coverage) in enumerate(top_anchor_kmers):
        if anchor_kmer not in anchor_groups:
            continue
            
        print(f"Processing anchor '{anchor_kmer}': {len(anchor_groups[anchor_kmer])} sequences", file=sys.stderr)
        
        anchor_processed_count = 0
        for seq_idx, anchor_pos_absolute in anchor_groups[anchor_kmer]:
            seq = sequences[seq_idx]
            header = headers[seq_idx]
            seq_str = str(seq).upper()
            
            # Calculate gaps to align anchor at global target position
            gaps_before = global_target_pos - anchor_pos_absolute
            
            # Limit gaps
            if gaps_before > max_gaps:
                print(f"  Skipping sequence {seq_idx}: needs {gaps_before} gaps (max {max_gaps})", file=sys.stderr)
                continue
            
            aligned_seq = '-' * gaps_before + seq_str
            
            # Calculate alignment score (higher score for better anchor rank)
            base_score = 1.0  # All anchored sequences get score 1.0
            anchor_boost = len(top_anchor_kmers) - anchor_rank
            score = base_score + anchor_boost
            
            aligned_sequences.append(aligned_seq)
            aligned_headers.append(header)
            alignment_scores.append(score)
            used_sequences.add(seq_idx)
            anchor_processed_count += 1
        
        print(f"  Added {anchor_processed_count} sequences with anchor '{anchor_kmer}'", file=sys.stderr)
    
    return aligned_headers, aligned_sequences, alignment_scores, used_sequences

def calculate_promoter_score(seq_str, anchor_kmer, anchor_pos):
    """Calculate alignment quality score based on anchor k-mer presence"""
    # Simple scoring based on having the anchor k-mer
    return 1.0 if anchor_kmer in seq_str else 0.0

def normalize_msa_length(aligned_sequences):
    """Make all sequences the same length by padding with gaps"""
    if not aligned_sequences:
        return aligned_sequences
    
    max_length = max(len(seq) for seq in aligned_sequences)
    normalized_sequences = []
    
    for seq in aligned_sequences:
        padded_seq = seq.ljust(max_length, '-')
        normalized_sequences.append(padded_seq)
    
    return normalized_sequences

def calculate_conservation_score(aligned_sequences, start, end):
    """Calculate conservation score for a region in the MSA"""
    if start >= end or end > len(aligned_sequences[0]):
        return 0.0
    
    total_score = 0.0
    total_positions = 0
    
    for pos in range(start, end):
        # Count nucleotides at this position (ignore gaps)
        nucleotides = {}
        valid_seqs = 0
        
        for seq in aligned_sequences:
            if pos < len(seq) and seq[pos] != '-':
                nucleotide = seq[pos].upper()
                nucleotides[nucleotide] = nucleotides.get(nucleotide, 0) + 1
                valid_seqs += 1
        
        if valid_seqs > 1:  # Need at least 2 sequences for conservation
            # Calculate conservation as the fraction of most common nucleotide
            max_count = max(nucleotides.values()) if nucleotides else 0
            conservation = max_count / valid_seqs if valid_seqs > 0 else 0.0
            total_score += conservation
            total_positions += 1
    
    return total_score / total_positions if total_positions > 0 else 0.0

def find_best_promoter_window(aligned_sequences, kmer_start, kmer_end, promoter_length):
    """Find the best conserved window of given length that contains the k-mer"""
    if not aligned_sequences:
        return kmer_start, kmer_end
    
    msa_length = len(aligned_sequences[0])
    kmer_length = kmer_end - kmer_start
    
    if promoter_length <= kmer_length:
        return kmer_start, kmer_end
    
    best_score = -1.0
    best_start = kmer_start
    best_end = kmer_end
    
    # Try all possible windows of promoter_length that contain the k-mer
    earliest_start = max(0, kmer_end - promoter_length)
    latest_start = min(msa_length - promoter_length, kmer_start)
    
    for window_start in range(earliest_start, latest_start + 1):
        window_end = window_start + promoter_length
        
        # Verify window contains the k-mer
        if window_start <= kmer_start and kmer_end <= window_end:
            score = calculate_conservation_score(aligned_sequences, window_start, window_end)
            
            if score > best_score:
                best_score = score
                best_start = window_start
                best_end = window_end
    
    print(f"Best promoter window: positions {best_start}-{best_end} (score: {best_score:.3f})", file=sys.stderr)
    return best_start, best_end

def annotate_promoter_regions(aligned_sequences, anchor_kmer, kmer_positions, promoter_length=20):
    """Annotate promoter regions in MSA with uppercase letters"""
    if not aligned_sequences:
        return aligned_sequences
    
    # Find consensus anchor position in the MSA
    consensus_anchor_pos = 0
    anchor_found = False
    
    # Look for anchor k-mer in first sequence to determine MSA coordinates
    first_seq = aligned_sequences[0].replace('-', '')
    if anchor_kmer in first_seq:
        # Find where anchor starts in ungapped sequence
        anchor_start_in_ungapped = first_seq.find(anchor_kmer)
        
        # Find corresponding position in gapped sequence
        ungapped_pos = 0
        for msa_pos, char in enumerate(aligned_sequences[0]):
            if char != '-':
                if ungapped_pos == anchor_start_in_ungapped:
                    consensus_anchor_pos = msa_pos
                    anchor_found = True
                    break
                ungapped_pos += 1
    
    if not anchor_found:
        print("Warning: Could not locate anchor in MSA, using default promoter region", file=sys.stderr)
        consensus_anchor_pos = 30  # Default fallback position
    
    # Calculate promoter region boundaries using sliding window approach
    kmer_length = len(anchor_kmer)
    kmer_start = consensus_anchor_pos
    kmer_end = consensus_anchor_pos + kmer_length
    
    if promoter_length == kmer_length:
        # If promoter length equals k-mer length, annotate exactly the k-mer positions
        promoter_start = kmer_start
        promoter_end = kmer_end
    else:
        # Use sliding window to find best conserved region containing the k-mer
        promoter_start, promoter_end = find_best_promoter_window(
            aligned_sequences, kmer_start, kmer_end, promoter_length)
    
    print(f"Annotating promoter region: MSA positions {promoter_start}-{promoter_end}", file=sys.stderr)
    
    # Annotate promoter regions
    annotated_sequences = []
    for seq in aligned_sequences:
        seq_list = list(seq.lower())
        
        # Convert promoter region to uppercase (gaps remain as gaps)
        for pos in range(promoter_start, promoter_end):
            if pos < len(seq_list):
                if seq_list[pos] != '-':
                    seq_list[pos] = seq_list[pos].upper()
        
        annotated_sequences.append(''.join(seq_list))
    
    return annotated_sequences

def add_rejected_sequences(sequences, headers, processed_indices, alignment_scores):
    """Add sequences that weren't included in MSA"""
    rejected_sequences = []
    rejected_headers = []
    
    for i, (seq, header) in enumerate(zip(sequences, headers)):
        if i not in processed_indices:
            rejected_sequences.append(str(seq).lower())
            rejected_headers.append(header)
    
    return rejected_headers, rejected_sequences

def main():
    parser = argparse.ArgumentParser(
        description='K-mer based promoter region detection and MSA',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default parameters
  python3 annotate_promoter_regions.py -i upstream.fa -o promoters_msa.fa
  
  # Custom k-mer size and promoter length
  python3 annotate_promoter_regions.py -i upstream.fa -k 8 -l 25 -o promoters_msa.fa
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input FASTA file with upstream sequences')
    parser.add_argument('-o', '--output', required=True,
                       help='Output MSA file with annotated promoter regions')
    parser.add_argument('-k', '--kmer-size', type=int, default=6,
                       help='K-mer size for anchor detection (default: 6)')
    parser.add_argument('-l', '--promoter-length', type=int, default=20,
                       help='Promoter region length (default: 20)')
    parser.add_argument('-n', '--max-gaps', type=int, default=40,
                       help='Maximum gaps at sequence ends (default: 40)')
    parser.add_argument('--promoter-start', type=int, default=-60,
                       help='Promoter region start relative to TSS (default: -60)')
    parser.add_argument('--promoter-end', type=int, default=-40,
                       help='Promoter region end relative to TSS (default: -40)')
    parser.add_argument('--min-sequences', type=int, default=3,
                       help='Minimum sequences required for anchor k-mer (default: 3)')
    parser.add_argument('-t', '--top-anchors', type=int, default=3,
                       help='Number of top anchor k-mers to use (default: 3)')
    
    args = parser.parse_args()
    
    # Load sequences
    print(f"Loading sequences from {args.input}...", file=sys.stderr)
    sequences = []
    headers = []
    
    try:
        with open(args.input, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                sequences.append(record.seq)
                headers.append(record.id)
        
        print(f"Loaded {len(sequences)} sequences", file=sys.stderr)
        
        if len(sequences) < args.min_sequences:
            print(f"Error: Need at least {args.min_sequences} sequences for analysis", file=sys.stderr)
            sys.exit(1)
        
        # Extract k-mers from promoter region
        print(f"Extracting {args.kmer_size}-mers from promoter region {args.promoter_start} to {args.promoter_end}...", file=sys.stderr)
        kmer_positions = extract_promoter_kmers(sequences, args.kmer_size, 
                                               args.promoter_start, args.promoter_end)
        
        if not kmer_positions:
            print("Error: No k-mers found in promoter regions", file=sys.stderr)
            sys.exit(1)
        
        # Find top t anchor k-mers
        top_anchor_kmers = find_best_anchor_kmers(
            kmer_positions, args.min_sequences, args.top_anchors)
        
        if not top_anchor_kmers:
            print(f"Error: No suitable anchor k-mers found (need k-mer in ≥{args.min_sequences} sequences)", file=sys.stderr)
            sys.exit(1)
        
        print(f"Found {len(top_anchor_kmers)} anchor k-mers:", file=sys.stderr)
        for i, (kmer, score, coverage) in enumerate(top_anchor_kmers):
            print(f"  #{i+1}: '{kmer}' (score: {score:.2f}, coverage: {coverage})", file=sys.stderr)
        
        # Align sequences using multiple anchors
        aligned_headers, aligned_sequences, alignment_scores, used_sequences = align_sequences_with_multiple_anchors(
            sequences, headers, top_anchor_kmers, kmer_positions, 
            args.promoter_length, args.max_gaps)
        
        if not aligned_sequences:
            print("Error: No sequences could be aligned with anchor k-mers", file=sys.stderr)
            sys.exit(1)
        
        # Normalize MSA length
        aligned_sequences = normalize_msa_length(aligned_sequences)
        
        # Annotate promoter regions (use first anchor for region detection)
        best_anchor_kmer = top_anchor_kmers[0][0]
        annotated_sequences = annotate_promoter_regions(
            aligned_sequences, best_anchor_kmer, kmer_positions, args.promoter_length)
        
        # Processed sequences are those that were assigned to anchors
        processed_indices = used_sequences
        
        # Add rejected sequences
        rejected_headers, rejected_sequences = add_rejected_sequences(
            sequences, headers, processed_indices, alignment_scores)
        
        # Include properly aligned sequences and rejected sequences at the end
        if annotated_sequences:
            max_length = max(len(seq) for seq in annotated_sequences)
            normalized_sequences = []
            for seq in annotated_sequences:
                normalized_sequences.append(seq.ljust(max_length, '-'))
        else:
            normalized_sequences = []
            max_length = 0
        
        # Add rejected sequences with end gaps to match MSA length
        rejected_sequences_padded = []
        for seq in rejected_sequences:
            # Convert to lowercase and pad with end gaps
            padded_seq = str(seq).lower().ljust(max_length, '-')
            rejected_sequences_padded.append(padded_seq)
        
        # Combine aligned and rejected sequences
        all_sequences = normalized_sequences + rejected_sequences_padded
        all_headers = list(aligned_headers) + list(rejected_headers)
        # Aligned sequences get their scores, rejected sequences get score 0
        all_scores = list(alignment_scores) + [0.0] * len(rejected_sequences)
        
        # Write output (no comments)
        print(f"Writing MSA to {args.output}...", file=sys.stderr)
        with open(args.output, 'w') as f:
            # Combine and sort all sequences by score (best first)
            sequence_data = list(zip(all_headers, all_sequences, all_scores))
            sequence_data.sort(key=lambda x: x[2], reverse=True)
            
            for header, seq, score in sequence_data:
                f.write(f">{header}\n{seq}\n")
        
        print(f"Done! MSA with {len(aligned_sequences)} aligned sequences + {len(rejected_sequences)} rejected (total: {len(all_sequences)})", file=sys.stderr)
        print(f"Promoter regions (length {args.promoter_length}) annotated in uppercase around {len(top_anchor_kmers)} anchor k-mers", file=sys.stderr)
        
    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()