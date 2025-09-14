#!/usr/bin/env python3
"""
Build Position Weight Matrices (PWMs) from aligned sequences.
Separates Pol II and Pol III sequences and builds PWMs with information vectors.
"""

import sys
import json
import argparse
import yaml
import numpy as np
from Bio import SeqIO
from collections import defaultdict
from pathlib import Path
import math

def load_polymerase_config(config_file: str = "config.yaml") -> dict:
    """Load RNA polymerase classification from config"""
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        return config.get('rna_polymerases', {})
    except Exception as e:
        print(f"Error loading config: {e}", file=sys.stderr)
        return {}

def parse_ncrna_from_header(header: str) -> str:
    """Extract ncRNA type from FASTA header"""
    # Expected format: >species-id|ncRNA_type|...
    parts = header.split('|')
    if len(parts) >= 2:
        return parts[1]
    return None

def extract_promoter_region(sequence: str, promoter_start: int = 40, promoter_end: int = 62) -> str:
    """
    Extract promoter region from aligned sequence.
    Default positions based on Drosophila analysis (columns 40-62).
    """
    # Remove gaps to check actual sequence length
    seq_no_gaps = sequence.replace('-', '')

    if len(seq_no_gaps) == 0:
        return None

    # Extract promoter region (keeping gaps for alignment)
    if len(sequence) > promoter_end:
        return sequence[promoter_start:promoter_end]
    else:
        return None

def build_pwm(sequences: list, pseudocount: float = 0.25) -> dict:
    """
    Build Position Weight Matrix from aligned sequences.
    Returns PWM with counts, frequencies, and information vector.
    """
    if not sequences:
        return None

    # Get length of alignment (should be same for all)
    length = len(sequences[0])

    # Initialize count matrix
    counts = {
        'A': [pseudocount] * length,
        'C': [pseudocount] * length,
        'G': [pseudocount] * length,
        'T': [pseudocount] * length
    }

    # Count nucleotides at each position
    for seq in sequences:
        for i, nt in enumerate(seq.upper()):
            if nt in 'ACGT':
                counts[nt][i] += 1

    # Convert counts to frequencies
    frequencies = {'A': [], 'C': [], 'G': [], 'T': []}
    for i in range(length):
        total = sum(counts[base][i] for base in 'ACGT')
        for base in 'ACGT':
            frequencies[base].append(counts[base][i] / total)

    # Calculate information vector
    information = []
    for i in range(length):
        info = 0
        for base in 'ACGT':
            f = frequencies[base][i]
            if f > 0:
                # Information content: f * log2(f/0.25)
                # Using log2(4*f) = log2(f/0.25)
                info += f * math.log2(4 * f)
        information.append(info)

    # Find core positions (5 most conserved consecutive positions)
    max_core_info = 0
    core_start = 0
    for i in range(length - 4):
        core_info = sum(information[i:i+5])
        if core_info > max_core_info:
            max_core_info = core_info
            core_start = i

    # Build PWM dictionary
    pwm = {
        'length': length,
        'sequences': len(sequences),
        'counts': counts,
        'frequencies': frequencies,
        'information': information,
        'core_start': core_start,
        'core_end': core_start + 5,
        'consensus': get_consensus(frequencies)
    }

    return pwm

def get_consensus(frequencies: dict) -> str:
    """Get consensus sequence from frequency matrix"""
    consensus = []
    length = len(frequencies['A'])

    for i in range(length):
        max_freq = 0
        max_base = 'N'
        for base in 'ACGT':
            if frequencies[base][i] > max_freq:
                max_freq = frequencies[base][i]
                max_base = base

        # Use IUPAC codes if frequency < 0.75
        if max_freq < 0.75:
            # Find second highest
            freqs = [(frequencies[b][i], b) for b in 'ACGT']
            freqs.sort(reverse=True)

            if freqs[1][0] > 0.25:  # Second base is significant
                bases = {freqs[0][1], freqs[1][1]}
                if bases == {'A', 'G'}:
                    max_base = 'R'
                elif bases == {'C', 'T'}:
                    max_base = 'Y'
                elif bases == {'A', 'T'}:
                    max_base = 'W'
                elif bases == {'G', 'C'}:
                    max_base = 'S'
                elif bases == {'G', 'T'}:
                    max_base = 'K'
                elif bases == {'A', 'C'}:
                    max_base = 'M'

        consensus.append(max_base)

    return ''.join(consensus)

def separate_sequences_by_polymerase(sequences: list, config_file: str) -> tuple:
    """Separate sequences into Pol II and Pol III based on ncRNA type"""
    polymerase_config = load_polymerase_config(config_file)

    pol2_types = set(polymerase_config.get('pol2', {}).get('types', []))
    pol3_types = set(polymerase_config.get('pol3', {}).get('types', []))

    pol2_sequences = []
    pol3_sequences = []

    for record in sequences:
        ncrna_type = parse_ncrna_from_header(record.id)

        if ncrna_type in pol2_types:
            pol2_sequences.append(str(record.seq))
        elif ncrna_type in pol3_types:
            pol3_sequences.append(str(record.seq))

    return pol2_sequences, pol3_sequences

def main():
    parser = argparse.ArgumentParser(
        description='Build PWMs from aligned sequences'
    )
    parser.add_argument('input', help='Input aligned FASTA file (annotated.fa)')
    parser.add_argument('-o', '--output-prefix', default='pwm',
                       help='Output prefix for PWM files (default: pwm)')
    parser.add_argument('-c', '--config', default='config.yaml',
                       help='Configuration file (default: config.yaml)')
    parser.add_argument('-s', '--promoter-start', type=int, default=40,
                       help='Start position of promoter in alignment (default: 40)')
    parser.add_argument('-e', '--promoter-end', type=int, default=62,
                       help='End position of promoter in alignment (default: 62)')
    parser.add_argument('-p', '--pseudocount', type=float, default=0.25,
                       help='Pseudocount for PWM (default: 0.25)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    # Read aligned sequences
    if args.verbose:
        print(f"Reading aligned sequences from {args.input}...", file=sys.stderr)

    sequences = list(SeqIO.parse(args.input, "fasta"))

    if not sequences:
        print("Error: No sequences found in input file", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print(f"Found {len(sequences)} sequences", file=sys.stderr)

    # Separate by polymerase
    pol2_seqs, pol3_seqs = separate_sequences_by_polymerase(sequences, args.config)

    if args.verbose:
        print(f"Pol II sequences: {len(pol2_seqs)}", file=sys.stderr)
        print(f"Pol III sequences: {len(pol3_seqs)}", file=sys.stderr)

    # Extract promoter regions
    pol2_promoters = []
    pol3_promoters = []

    for seq in pol2_seqs:
        promoter = extract_promoter_region(seq, args.promoter_start, args.promoter_end)
        if promoter:
            pol2_promoters.append(promoter)

    for seq in pol3_seqs:
        promoter = extract_promoter_region(seq, args.promoter_start, args.promoter_end)
        if promoter:
            pol3_promoters.append(promoter)

    if args.verbose:
        print(f"Pol II promoters extracted: {len(pol2_promoters)}", file=sys.stderr)
        print(f"Pol III promoters extracted: {len(pol3_promoters)}", file=sys.stderr)

    # Build PWMs
    if pol2_promoters:
        pol2_pwm = build_pwm(pol2_promoters, args.pseudocount)
        if pol2_pwm:
            output_file = f"{args.output_prefix}_pol2.json"
            with open(output_file, 'w') as f:
                json.dump(pol2_pwm, f, indent=2)
            print(f"Pol II PWM written to {output_file}", file=sys.stderr)
            print(f"  Consensus: {pol2_pwm['consensus']}", file=sys.stderr)
            print(f"  Core positions: {pol2_pwm['core_start']}-{pol2_pwm['core_end']}", file=sys.stderr)

    if pol3_promoters:
        pol3_pwm = build_pwm(pol3_promoters, args.pseudocount)
        if pol3_pwm:
            output_file = f"{args.output_prefix}_pol3.json"
            with open(output_file, 'w') as f:
                json.dump(pol3_pwm, f, indent=2)
            print(f"Pol III PWM written to {output_file}", file=sys.stderr)
            print(f"  Consensus: {pol3_pwm['consensus']}", file=sys.stderr)
            print(f"  Core positions: {pol3_pwm['core_start']}-{pol3_pwm['core_end']}", file=sys.stderr)

if __name__ == "__main__":
    main()