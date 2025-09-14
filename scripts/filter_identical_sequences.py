#!/usr/bin/env python3
"""
Filter out nearly identical sequences from FASTA input
Reads FASTA from stdin and removes sequences that are above a similarity threshold
"""

import sys
import argparse
from typing import Dict, List, Set
from collections import defaultdict

def parse_fasta(file_handle):
    """Parse FASTA sequences from file handle"""
    sequences = {}
    current_id = None
    current_seq = []
    
    for line in file_handle:
        line = line.strip()
        if line.startswith('>'):
            # Save previous sequence
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
            # Start new sequence
            current_id = line[1:]  # Remove '>' character
            current_seq = []
        else:
            current_seq.append(line.upper())
    
    # Save last sequence
    if current_id is not None:
        sequences[current_id] = ''.join(current_seq)
    
    return sequences

def calculate_similarity(seq1: str, seq2: str) -> float:
    """Calculate percentage similarity between two sequences"""
    if len(seq1) != len(seq2):
        return 0.0
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / len(seq1)) * 100.0

def find_redundant_sequences(sequences: Dict[str, str], similarity_threshold: float) -> Set[str]:
    """Find sequences that are above similarity threshold with other sequences"""
    seq_ids = list(sequences.keys())
    redundant = set()
    
    # Group sequences by length for efficiency
    length_groups = defaultdict(list)
    for seq_id, seq in sequences.items():
        length_groups[len(seq)].append(seq_id)
    
    # Only compare sequences of the same length
    for length, ids_in_group in length_groups.items():
        for i in range(len(ids_in_group)):
            if ids_in_group[i] in redundant:
                continue
                
            seq1 = sequences[ids_in_group[i]]
            
            for j in range(i + 1, len(ids_in_group)):
                if ids_in_group[j] in redundant:
                    continue
                    
                seq2 = sequences[ids_in_group[j]]
                similarity = calculate_similarity(seq1, seq2)
                
                if similarity >= similarity_threshold:
                    # Keep the first sequence, mark the second as redundant
                    redundant.add(ids_in_group[j])
    
    return redundant

def filter_sequences(sequences: Dict[str, str], similarity_threshold: float) -> Dict[str, str]:
    """Filter out redundant sequences based on similarity threshold"""
    redundant = find_redundant_sequences(sequences, similarity_threshold)
    
    # Return sequences that are not redundant
    filtered = {seq_id: seq for seq_id, seq in sequences.items() if seq_id not in redundant}
    
    return filtered, redundant

def main():
    parser = argparse.ArgumentParser(
        description='Filter out nearly identical sequences from FASTA input',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script removes sequences that are above a specified similarity threshold.
When two sequences are similar, the first one encountered is kept.

Examples:
  # Filter out sequences that are 95% or more identical
  cat sequences.fa | python3 filter_identical_sequences.py -s 95 > filtered.fa
  
  # Filter out sequences that are 99% or more identical (very strict)
  cat sequences.fa | python3 filter_identical_sequences.py -s 99 > filtered.fa
  
  # In a pipeline
  cat upstream.fa | python3 filter_identical_sequences.py -s 99 | other_processing.py

Note: Only sequences of the same length are compared for efficiency.
      Sequences are compared character by character for exact similarity calculation.
        """
    )
    
    parser.add_argument('-s', '--similarity', type=float, default=99.0,
                       help='Similarity threshold percentage (default: 99.0)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Print filtering statistics to stderr')
    
    args = parser.parse_args()
    
    if not (0 <= args.similarity <= 100):
        print("Error: Similarity threshold must be between 0 and 100", file=sys.stderr)
        sys.exit(1)
    
    # Read sequences from stdin
    sequences = parse_fasta(sys.stdin)
    
    if len(sequences) == 0:
        print("Warning: No sequences found in input", file=sys.stderr)
        sys.exit(0)
    
    # Filter sequences
    filtered_sequences, redundant_ids = filter_sequences(sequences, args.similarity)
    
    # Always report statistics to stderr (not just in verbose mode)
    removed_count = len(redundant_ids)
    kept_count = len(filtered_sequences)
    print(f"Read {len(sequences)} sequences", file=sys.stderr)
    print(f"Filtered {removed_count} sequences above {args.similarity}% similarity", file=sys.stderr)
    print(f"Kept {kept_count} unique sequences", file=sys.stderr)
    if len(sequences) > 0:
        print(f"Reduction: {removed_count/len(sequences)*100:.1f}%", file=sys.stderr)
    
    # Output filtered sequences to stdout
    for seq_id, sequence in filtered_sequences.items():
        print(f">{seq_id}")
        # Print sequence in 80-character lines
        for i in range(0, len(sequence), 80):
            print(sequence[i:i+80])

if __name__ == "__main__":
    main()