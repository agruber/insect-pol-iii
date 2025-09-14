#!/usr/bin/env python3
"""
Filter out low complexity sequences from FASTA input
Removes sequences with excessive homopolymers, dinucleotide or trinucleotide repeats
"""

import sys
import argparse
from collections import Counter
from typing import Dict, Tuple

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

def check_homopolymer_runs(sequence: str, threshold: float = 0.5) -> bool:
    """
    Check if any single nucleotide comprises more than threshold of the sequence
    
    Args:
        sequence: DNA sequence
        threshold: Maximum fraction of sequence for any single nucleotide
        
    Returns:
        True if sequence fails (is low complexity), False if passes
    """
    if not sequence:
        return False
    
    nucleotide_counts = Counter(sequence)
    seq_length = len(sequence)
    
    for nucleotide, count in nucleotide_counts.items():
        if count / seq_length > threshold:
            return True  # Low complexity detected
    
    return False

def check_dinucleotide_repeats(sequence: str, threshold: float = 0.75) -> bool:
    """
    Check if any dinucleotide pattern comprises more than threshold of the sequence
    
    Args:
        sequence: DNA sequence
        threshold: Maximum fraction of sequence for any dinucleotide pattern
        
    Returns:
        True if sequence fails (is low complexity), False if passes
    """
    if len(sequence) < 2:
        return False
    
    seq_length = len(sequence)
    
    # Check all possible dinucleotides
    dinucleotides = ['AA', 'AT', 'AC', 'AG',
                     'TA', 'TT', 'TC', 'TG',
                     'CA', 'CT', 'CC', 'CG',
                     'GA', 'GT', 'GC', 'GG']
    
    for dinuc in dinucleotides:
        # Count occurrences of the dinucleotide pattern
        count = 0
        for i in range(0, len(sequence) - 1, 2):
            if sequence[i:i+2] == dinuc:
                count += 1
        
        # Check if this dinucleotide dominates the sequence
        if (count * 2) / seq_length > threshold:
            return True  # Low complexity detected
    
    # Also check shifted dinucleotides (starting at position 1)
    for dinuc in dinucleotides:
        count = 0
        for i in range(1, len(sequence) - 1, 2):
            if sequence[i:i+2] == dinuc:
                count += 1
        
        if (count * 2) / seq_length > threshold:
            return True  # Low complexity detected
    
    return False

def check_trinucleotide_repeats(sequence: str, threshold: float = 0.75) -> bool:
    """
    Check if any trinucleotide pattern comprises more than threshold of the sequence
    
    Args:
        sequence: DNA sequence
        threshold: Maximum fraction of sequence for any trinucleotide pattern
        
    Returns:
        True if sequence fails (is low complexity), False if passes
    """
    if len(sequence) < 3:
        return False
    
    seq_length = len(sequence)
    
    # Count all trinucleotides in the sequence
    trinuc_counts = {}
    
    # Check all three reading frames
    for frame in range(3):
        for i in range(frame, len(sequence) - 2, 3):
            trinuc = sequence[i:i+3]
            if len(trinuc) == 3:  # Ensure we have a full trinucleotide
                if trinuc not in trinuc_counts:
                    trinuc_counts[trinuc] = 0
                trinuc_counts[trinuc] += 1
    
    # Check if any trinucleotide pattern is overrepresented
    for trinuc, count in trinuc_counts.items():
        if (count * 3) / seq_length > threshold:
            return True  # Low complexity detected
    
    return False

def is_low_complexity(sequence: str, homopolymer_threshold: float = 0.5,
                     dinuc_threshold: float = 0.75, trinuc_threshold: float = 0.75) -> Tuple[bool, str]:
    """
    Check if a sequence is low complexity based on multiple criteria
    
    Returns:
        Tuple of (is_low_complexity, reason)
    """
    # Check homopolymer runs
    if check_homopolymer_runs(sequence, homopolymer_threshold):
        return True, "homopolymer"
    
    # Check dinucleotide repeats
    if check_dinucleotide_repeats(sequence, dinuc_threshold):
        return True, "dinucleotide"
    
    # Check trinucleotide repeats
    if check_trinucleotide_repeats(sequence, trinuc_threshold):
        return True, "trinucleotide"
    
    return False, ""

def filter_sequences(sequences: Dict[str, str], homopolymer_threshold: float = 0.5,
                    dinuc_threshold: float = 0.75, trinuc_threshold: float = 0.75) -> Dict[str, str]:
    """Filter out low complexity sequences"""
    filtered = {}
    removed_counts = {'homopolymer': 0, 'dinucleotide': 0, 'trinucleotide': 0}
    
    for seq_id, sequence in sequences.items():
        is_low, reason = is_low_complexity(sequence, homopolymer_threshold, 
                                          dinuc_threshold, trinuc_threshold)
        if not is_low:
            filtered[seq_id] = sequence
        else:
            removed_counts[reason] += 1
    
    return filtered, removed_counts

def main():
    parser = argparse.ArgumentParser(
        description='Filter out low complexity sequences from FASTA input',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script removes sequences with low complexity patterns:

1. Homopolymer runs: Any single nucleotide >50% of sequence
2. Dinucleotide repeats: Any dinucleotide pattern >75% of sequence  
3. Trinucleotide repeats: Any trinucleotide pattern >75% of sequence

Examples:
  # Filter with default thresholds
  cat sequences.fa | python3 filter_low_complexity.py > filtered.fa
  
  # With custom homopolymer threshold
  cat sequences.fa | python3 filter_low_complexity.py --homopolymer 0.6 > filtered.fa
  
  # With custom dinucleotide threshold
  cat sequences.fa | python3 filter_low_complexity.py --dinucleotide 0.8 > filtered.fa

Note: All thresholds are expressed as fractions (0.5 = 50%, 0.75 = 75%)
        """
    )
    
    parser.add_argument('--homopolymer', type=float, default=0.5,
                       help='Maximum fraction for any single nucleotide (default: 0.5)')
    parser.add_argument('--dinucleotide', type=float, default=0.75,
                       help='Maximum fraction for any dinucleotide pattern (default: 0.75)')
    parser.add_argument('--trinucleotide', type=float, default=0.75,
                       help='Maximum fraction for any trinucleotide pattern (default: 0.75)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Print detailed filtering statistics to stderr')
    
    args = parser.parse_args()
    
    # Validate thresholds
    for threshold, name in [(args.homopolymer, 'homopolymer'), 
                           (args.dinucleotide, 'dinucleotide'),
                           (args.trinucleotide, 'trinucleotide')]:
        if not (0 < threshold <= 1):
            print(f"Error: {name} threshold must be between 0 and 1", file=sys.stderr)
            sys.exit(1)
    
    # Read sequences from stdin
    sequences = parse_fasta(sys.stdin)
    
    if len(sequences) == 0:
        print("Warning: No sequences found in input", file=sys.stderr)
        sys.exit(0)
    
    # Filter sequences
    filtered_sequences, removed_counts = filter_sequences(
        sequences, args.homopolymer, args.dinucleotide, args.trinucleotide
    )
    
    # Report statistics to stderr
    total_removed = sum(removed_counts.values())
    print(f"Read {len(sequences)} sequences", file=sys.stderr)
    print(f"Filtered {total_removed} low complexity sequences", file=sys.stderr)
    if total_removed > 0:
        print(f"  Homopolymer runs: {removed_counts['homopolymer']}", file=sys.stderr)
        print(f"  Dinucleotide repeats: {removed_counts['dinucleotide']}", file=sys.stderr)
        print(f"  Trinucleotide repeats: {removed_counts['trinucleotide']}", file=sys.stderr)
    print(f"Kept {len(filtered_sequences)} sequences", file=sys.stderr)
    if len(sequences) > 0:
        print(f"Reduction: {total_removed/len(sequences)*100:.1f}%", file=sys.stderr)
    
    # Output filtered sequences to stdout
    for seq_id, sequence in filtered_sequences.items():
        print(f">{seq_id}")
        # Print sequence in 80-character lines
        for i in range(0, len(sequence), 80):
            print(sequence[i:i+80])

if __name__ == "__main__":
    main()