#!/usr/bin/env python3
"""
Implementation of the MATCH™ algorithm for transcription factor binding site prediction.
Based on Kel et al. 2003.
"""

import math
import json
from typing import Dict, List, Tuple
from collections import defaultdict

class MATCHAlgorithm:
    """
    MATCH™ algorithm implementation with information vector weighting.
    """

    def __init__(self, pwm: dict, cutoff_strategy: str = 'minSum'):
        """
        Initialize MATCH algorithm with a PWM.

        Args:
            pwm: Position Weight Matrix dictionary with frequencies and information vector
            cutoff_strategy: One of 'minFN', 'minFP', or 'minSum'
        """
        self.pwm = pwm
        self.length = pwm['length']
        self.frequencies = pwm['frequencies']
        self.information = pwm['information']
        self.core_start = pwm['core_start']
        self.core_end = pwm['core_end']

        # Pre-calculate min and max scores for normalization
        self._calculate_score_ranges()

        # Set cutoff thresholds based on strategy
        self._set_cutoffs(cutoff_strategy)

        # Build hash table for optimization
        self.pentanucleotide_scores = {}

    def _calculate_score_ranges(self):
        """Calculate minimum and maximum possible scores for the PWM"""
        self.min_scores = []
        self.max_scores = []

        for i in range(self.length):
            # Find min and max frequency at this position
            freqs = [self.frequencies[base][i] for base in 'ACGT']
            min_freq = min(freqs)
            max_freq = max(freqs)

            # Weight by information content
            self.min_scores.append(self.information[i] * math.log2(4 * min_freq) if min_freq > 0 else -10)
            self.max_scores.append(self.information[i] * math.log2(4 * max_freq) if max_freq > 0 else 0)

        self.total_min = sum(self.min_scores)
        self.total_max = sum(self.max_scores)

        # Same for core region
        self.core_min = sum(self.min_scores[self.core_start:self.core_end])
        self.core_max = sum(self.max_scores[self.core_start:self.core_end])

    def _set_cutoffs(self, strategy: str):
        """Set CSS and MSS cutoff thresholds based on strategy"""
        if strategy == 'minFN':
            # Minimize false negatives - more sensitive
            self.css_cutoff = 0.60  # Accept 60% core similarity
            self.mss_cutoff = 0.70  # Accept 70% matrix similarity
        elif strategy == 'minFP':
            # Minimize false positives - more specific
            self.css_cutoff = 0.85  # Require 85% core similarity
            self.mss_cutoff = 0.90  # Require 90% matrix similarity
        else:  # minSum - balanced
            self.css_cutoff = 0.75  # Require 75% core similarity
            self.mss_cutoff = 0.80  # Require 80% matrix similarity

    def calculate_score(self, sequence: str, start: int, end: int) -> float:
        """
        Calculate information-weighted score for a sequence region.

        Args:
            sequence: DNA sequence
            start: Start position in PWM
            end: End position in PWM

        Returns:
            Raw score (not normalized)
        """
        score = 0
        seq_upper = sequence.upper()

        for i in range(start, min(end, len(seq_upper))):
            if i >= len(seq_upper):
                break

            nt = seq_upper[i]
            if nt in 'ACGT':
                pwm_pos = i - start
                if pwm_pos < self.length:
                    freq = self.frequencies[nt][pwm_pos]
                    if freq > 0:
                        # Information-weighted score
                        score += self.information[pwm_pos] * math.log2(4 * freq)
                    else:
                        score += self.min_scores[pwm_pos]

        return score

    def calculate_css(self, pentamer: str) -> float:
        """
        Calculate Core Similarity Score for a pentanucleotide.

        Args:
            pentamer: 5-mer sequence

        Returns:
            CSS normalized to [0, 1]
        """
        if len(pentamer) != 5:
            return 0

        score = 0
        pentamer_upper = pentamer.upper()

        for i in range(5):
            nt = pentamer_upper[i]
            if nt in 'ACGT':
                pwm_pos = self.core_start + i
                if pwm_pos < self.length:
                    freq = self.frequencies[nt][pwm_pos]
                    if freq > 0:
                        score += self.information[pwm_pos] * math.log2(4 * freq)
                    else:
                        score += self.min_scores[pwm_pos]

        # Normalize to [0, 1]
        if self.core_max - self.core_min > 0:
            css = (score - self.core_min) / (self.core_max - self.core_min)
        else:
            css = 0

        return max(0, min(1, css))  # Clamp to [0, 1]

    def calculate_mss(self, sequence: str) -> float:
        """
        Calculate Matrix Similarity Score for a sequence.

        Args:
            sequence: DNA sequence of PWM length

        Returns:
            MSS normalized to [0, 1]
        """
        if len(sequence) < self.length:
            return 0

        score = self.calculate_score(sequence, 0, self.length)

        # Normalize to [0, 1]
        if self.total_max - self.total_min > 0:
            mss = (score - self.total_min) / (self.total_max - self.total_min)
        else:
            mss = 0

        return max(0, min(1, mss))  # Clamp to [0, 1]

    def build_pentanucleotide_hash(self, sequence: str) -> Dict[str, List[int]]:
        """
        Build hash table of pentanucleotides in the sequence with their positions.

        Args:
            sequence: DNA sequence to scan

        Returns:
            Dictionary mapping 5-mers to list of positions
        """
        hash_table = defaultdict(list)
        seq_upper = sequence.upper()

        for i in range(len(seq_upper) - 4):
            pentamer = seq_upper[i:i+5]
            # Only process valid DNA pentamers
            if all(nt in 'ACGT' for nt in pentamer):
                hash_table[pentamer].append(i)

        return hash_table

    def scan_sequence(self, sequence: str, both_strands: bool = True) -> List[Dict]:
        """
        Scan a DNA sequence for matches to the PWM using MATCH algorithm.

        Args:
            sequence: DNA sequence to scan
            both_strands: Whether to scan both strands

        Returns:
            List of matches with positions and scores
        """
        matches = []

        # Build pentanucleotide hash table
        penta_hash = self.build_pentanucleotide_hash(sequence)

        # Pre-calculate CSS for all unique pentanucleotides
        penta_css = {}
        for pentamer in penta_hash.keys():
            # Check if pentamer aligns with core region
            css = self.calculate_css(pentamer)
            if css >= self.css_cutoff:
                penta_css[pentamer] = css

        # Scan using pentanucleotides that pass CSS threshold
        for pentamer, css in penta_css.items():
            for pos in penta_hash[pentamer]:
                # Determine where this pentamer aligns with the PWM core
                # Try different alignments of the pentamer with the core
                for core_offset in range(5):
                    pwm_start = pos - (self.core_start + core_offset)

                    if pwm_start >= 0 and pwm_start + self.length <= len(sequence):
                        # Extract full sequence for PWM
                        test_seq = sequence[pwm_start:pwm_start + self.length]

                        # Calculate MSS for full match
                        mss = self.calculate_mss(test_seq)

                        if mss >= self.mss_cutoff:
                            matches.append({
                                'start': pwm_start,
                                'end': pwm_start + self.length,
                                'strand': '+',
                                'sequence': test_seq,
                                'css': css,
                                'mss': mss,
                                'score': css * 0.3 + mss * 0.7  # Combined score
                            })

        # Scan reverse complement if requested
        if both_strands:
            rev_comp = self.reverse_complement(sequence)
            rev_matches = self.scan_sequence(rev_comp, both_strands=False)

            # Adjust positions for reverse strand
            for match in rev_matches:
                match['strand'] = '-'
                # Convert positions back to original strand coordinates
                original_start = len(sequence) - match['end']
                original_end = len(sequence) - match['start']
                match['start'] = original_start
                match['end'] = original_end

            matches.extend(rev_matches)

        # Remove duplicates and sort by score
        unique_matches = self._remove_overlapping_matches(matches)

        return sorted(unique_matches, key=lambda x: x['score'], reverse=True)

    def _remove_overlapping_matches(self, matches: List[Dict]) -> List[Dict]:
        """Remove overlapping matches, keeping the highest scoring ones"""
        if not matches:
            return []

        # Sort by score (descending)
        sorted_matches = sorted(matches, key=lambda x: x['score'], reverse=True)

        kept = []
        for match in sorted_matches:
            # Check if this match overlaps with any kept match
            overlaps = False
            for kept_match in kept:
                if match['strand'] == kept_match['strand']:
                    if (match['start'] < kept_match['end'] and
                        match['end'] > kept_match['start']):
                        overlaps = True
                        break

            if not overlaps:
                kept.append(match)

        return kept

    @staticmethod
    def reverse_complement(sequence: str) -> str:
        """Get reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                     'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                     'N': 'N', 'n': 'n'}
        return ''.join(complement.get(base, base) for base in sequence[::-1])

def load_pwm(pwm_file: str) -> dict:
    """Load PWM from JSON file"""
    with open(pwm_file, 'r') as f:
        return json.load(f)

def main():
    """Test the MATCH algorithm"""
    import argparse

    parser = argparse.ArgumentParser(description='Test MATCH algorithm')
    parser.add_argument('pwm_file', help='PWM JSON file')
    parser.add_argument('sequence', help='DNA sequence to scan')
    parser.add_argument('-c', '--cutoff', default='minSum',
                       choices=['minFN', 'minFP', 'minSum'],
                       help='Cutoff strategy (default: minSum)')

    args = parser.parse_args()

    # Load PWM
    pwm = load_pwm(args.pwm_file)

    # Initialize MATCH algorithm
    matcher = MATCHAlgorithm(pwm, args.cutoff)

    # Scan sequence
    matches = matcher.scan_sequence(args.sequence)

    # Print results
    print(f"Found {len(matches)} matches:")
    for i, match in enumerate(matches, 1):
        print(f"\nMatch {i}:")
        print(f"  Position: {match['start']}-{match['end']} ({match['strand']})")
        print(f"  Sequence: {match['sequence']}")
        print(f"  CSS: {match['css']:.3f}")
        print(f"  MSS: {match['mss']:.3f}")
        print(f"  Score: {match['score']:.3f}")

if __name__ == "__main__":
    main()