#!/usr/bin/env python3
"""
Analyze promoter consensus sequences from aligned promoter data.
Generates IUPAC consensus sequences and identifies Pol II/III discriminating positions.
"""

import sys
import os
import argparse
import yaml
from typing import Dict, List, Tuple, Set
from collections import Counter, defaultdict
from Bio import SeqIO
import numpy as np

def load_rna_config(config_file: str = "config.yaml") -> Dict:
    """Load RNA polymerase classification from config file"""
    config_path = os.path.join(os.path.dirname(__file__), '..', config_file)
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        return config.get('rna_polymerases', {})
    except Exception as e:
        print(f"Warning: Could not load config file: {e}", file=sys.stderr)
        return {}

def classify_sequence(seq_id: str, rna_config: Dict) -> str:
    """Classify sequence as pol2, pol3, or unknown based on RNA type"""
    # Extract RNA type from sequence ID (assuming format: species|type|...)
    parts = seq_id.split('|')
    if len(parts) < 2:
        return 'unknown'

    rna_type = parts[1]

    # Check pol2 types
    pol2_types = rna_config.get('pol2', {}).get('types', [])
    if rna_type in pol2_types:
        return 'pol2'

    # Check pol3 types
    pol3_types = rna_config.get('pol3', {}).get('types', [])
    if rna_type in pol3_types:
        return 'pol3'

    return 'unknown'

def extract_promoter_regions(sequences: List[str], seq_ids: List[str]) -> Tuple[List[str], List[int], List[int]]:
    """
    Extract uppercase promoter regions from aligned sequences.

    Returns:
        - promoter_sequences: List of promoter sequences (uppercase parts only)
        - start_positions: Start position of promoter in each sequence
        - end_positions: End position of promoter in each sequence
    """
    promoter_sequences = []
    start_positions = []
    end_positions = []

    for seq in sequences:
        # Find start and end of uppercase region
        start = -1
        end = -1

        for i, char in enumerate(seq):
            if char.isupper() and char != '-':
                if start == -1:
                    start = i
                end = i

        if start != -1:
            # Extract the promoter region (including gaps within it)
            promoter_region = seq[start:end+1]
            # Remove gaps from promoter sequence for consensus calculation
            promoter_clean = promoter_region.replace('-', '')
            promoter_sequences.append(promoter_clean)
            start_positions.append(start)
            end_positions.append(end)
        else:
            # No uppercase region found
            promoter_sequences.append("")
            start_positions.append(-1)
            end_positions.append(-1)

    return promoter_sequences, start_positions, end_positions

def calculate_nucleotide_frequencies(sequences: List[str]) -> Dict[int, Dict[str, float]]:
    """Calculate nucleotide frequencies at each position"""
    if not sequences or not sequences[0]:
        return {}

    # Find the maximum length among all sequences
    max_len = max(len(seq) for seq in sequences if seq)

    frequencies = {}

    for pos in range(max_len):
        pos_counts = Counter()
        total = 0

        for seq in sequences:
            if seq and pos < len(seq):
                base = seq[pos].upper()
                if base in 'ACGT':
                    pos_counts[base] += 1
                    total += 1

        if total > 0:
            frequencies[pos] = {
                base: count / total
                for base, count in pos_counts.items()
            }

    return frequencies

def frequency_to_iupac(freq_dict: Dict[str, float], threshold: float = 0.5) -> str:
    """Convert nucleotide frequencies to IUPAC code"""
    if not freq_dict:
        return 'N'

    # Sort by frequency (highest first)
    sorted_bases = sorted(freq_dict.items(), key=lambda x: x[1], reverse=True)

    # Determine which bases to include (above threshold)
    included_bases = set()
    for base, freq in sorted_bases:
        if freq >= threshold:
            included_bases.add(base)

    # If no base meets threshold, take the most frequent one
    if not included_bases:
        included_bases = {sorted_bases[0][0]}

    # Convert to IUPAC code
    bases_str = ''.join(sorted(included_bases))

    iupac_map = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K',
        'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B',
        'ACGT': 'N'
    }

    return iupac_map.get(bases_str, 'N')

def generate_consensus_sequence(sequences: List[str], threshold: float = 0.5) -> str:
    """Generate IUPAC consensus sequence from a list of sequences"""
    frequencies = calculate_nucleotide_frequencies(sequences)
    consensus = ""

    max_pos = max(frequencies.keys()) if frequencies else -1

    for pos in range(max_pos + 1):
        if pos in frequencies:
            iupac_char = frequency_to_iupac(frequencies[pos], threshold)
            consensus += iupac_char
        else:
            consensus += 'N'

    return consensus

def generate_position_wise_consensus(sequences: List[str], discriminatory_positions: set, threshold: float = 0.2) -> str:
    """
    Generate consensus sequence with position-specific logic:
    - Non-discriminatory positions: Use combined frequencies from all sequences
    - Discriminatory positions: Will be marked as '|' in display
    """
    frequencies = calculate_nucleotide_frequencies(sequences)
    consensus = ""

    max_pos = max(frequencies.keys()) if frequencies else -1

    for pos in range(max_pos + 1):
        if pos in frequencies:
            if pos in discriminatory_positions:
                # For discriminatory positions, we'll replace with '|' in display
                # But store the combined consensus here for non-display purposes
                iupac_char = frequency_to_iupac(frequencies[pos], threshold)
                consensus += iupac_char
            else:
                # Non-discriminatory positions: use combined frequencies
                iupac_char = frequency_to_iupac(frequencies[pos], threshold)
                consensus += iupac_char
        else:
            consensus += 'N'

    return consensus

def iupac_codes_overlap(code1: str, code2: str) -> bool:
    """Check if two IUPAC codes have overlapping nucleotides"""
    iupac_to_bases = {
        'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
        'M': {'A', 'C'}, 'R': {'A', 'G'}, 'W': {'A', 'T'},
        'S': {'C', 'G'}, 'Y': {'C', 'T'}, 'K': {'G', 'T'},
        'V': {'A', 'C', 'G'}, 'H': {'A', 'C', 'T'},
        'D': {'A', 'G', 'T'}, 'B': {'C', 'G', 'T'},
        'N': {'A', 'C', 'G', 'T'}
    }

    bases1 = iupac_to_bases.get(code1, set())
    bases2 = iupac_to_bases.get(code2, set())

    # Check if there's any overlap
    return bool(bases1 & bases2)

def calculate_jensen_shannon_divergence(p: Dict[str, float], q: Dict[str, float]) -> float:
    """
    Calculate Jensen-Shannon divergence between two distributions.
    JS divergence is symmetric and bounded to [0, 1] when using log base 2.
    """
    epsilon = 1e-10  # Small value to avoid log(0)

    # Create the average distribution
    all_bases = set(p.keys()) | set(q.keys())
    m = {}
    for base in all_bases:
        p_val = p.get(base, epsilon)
        q_val = q.get(base, epsilon)
        m[base] = (p_val + q_val) / 2.0

    # Calculate KL divergences
    kl_p_m = 0.0
    kl_q_m = 0.0

    for base in all_bases:
        p_val = p.get(base, epsilon)
        q_val = q.get(base, epsilon)
        m_val = m[base]

        if p_val > epsilon:
            kl_p_m += p_val * np.log2(p_val / m_val)  # Using log2 for bounded [0,1]
        if q_val > epsilon:
            kl_q_m += q_val * np.log2(q_val / m_val)

    # JS divergence is the average of the two KL divergences
    js_div = (kl_p_m + kl_q_m) / 2.0

    # JS divergence with log2 is bounded to [0, 1]
    return min(max(js_div, 0.0), 1.0)  # Ensure it's in [0, 1]

def calculate_position_discrimination(pol2_freqs: Dict[int, Dict[str, float]],
                                    pol3_freqs: Dict[int, Dict[str, float]]) -> List[Tuple[int, float, str, str]]:
    """
    Calculate discriminating power at each position using Jensen-Shannon divergence.
    JS divergence is normalized to [0, 1].

    Returns:
        List of (position, discrimination_score, pol2_consensus, pol3_consensus)
    """
    discriminating_positions = []

    # Get all positions that exist in both groups
    common_positions = set(pol2_freqs.keys()) & set(pol3_freqs.keys())

    for pos in sorted(common_positions):
        pol2_freq = pol2_freqs[pos]
        pol3_freq = pol3_freqs[pos]

        # Calculate JS divergence (automatically normalized to [0, 1])
        discrimination_score = calculate_jensen_shannon_divergence(pol2_freq, pol3_freq)

        # Generate consensus for each group at this position
        pol2_consensus = frequency_to_iupac(pol2_freq, 0.2)
        pol3_consensus = frequency_to_iupac(pol3_freq, 0.2)

        # Include all positions with significant divergence
        if discrimination_score > 0.5:  # Threshold for high divergence (50% of maximum)
            discriminating_positions.append((pos, discrimination_score, pol2_consensus, pol3_consensus))

    # Sort by discrimination score (highest first)
    discriminating_positions.sort(key=lambda x: x[1], reverse=True)

    return discriminating_positions

def format_discriminating_positions(disc_positions: List[Tuple[int, float, str, str]],
                                  max_positions: int = 10) -> str:
    """Format discriminating positions for display"""
    if not disc_positions:
        return "No discriminating positions found"

    # Sort by position number for readability
    sorted_positions = sorted(disc_positions[:max_positions], key=lambda x: x[0])

    formatted_parts = []
    for pos, score, pol2_base, pol3_base in sorted_positions:
        formatted_parts.append(f"pos{pos+1}({pol2_base}/{pol3_base}:JS={score:.3f})")

    return ", ".join(formatted_parts)

def create_compact_consensus_display(overall_consensus: str, pol2_consensus: str, pol3_consensus: str,
                                   disc_positions: List[Tuple[int, float, str, str]]) -> str:
    """Create compact consensus display with polymerase-specific nucleotides marked"""

    if not pol2_consensus or not pol3_consensus:
        return f"Consensus: {overall_consensus}"

    # Use all positions with JS > 0.5 (they're already filtered in calculate_position_discrimination)
    display_positions = []
    warnings = []

    for pos, score, pol2_base, pol3_base in disc_positions:
        display_positions.append((pos, pol2_base, pol3_base))
        # Check for IUPAC overlap and warn if present
        if iupac_codes_overlap(pol2_base, pol3_base):
            # Get the nucleotide explanations
            iupac_to_bases = {
                'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
                'M': 'A/C', 'R': 'A/G', 'W': 'A/T', 'S': 'C/G', 'Y': 'C/T', 'K': 'G/T',
                'V': 'A/C/G', 'H': 'A/C/T', 'D': 'A/G/T', 'B': 'C/G/T',
                'N': 'A/C/G/T'
            }
            pol2_exp = iupac_to_bases.get(pol2_base, pol2_base)
            pol3_exp = iupac_to_bases.get(pol3_base, pol3_base)
            warnings.append(f"Warning: Position {pos+1} has overlapping IUPAC codes ({pol2_base}={pol2_exp} vs {pol3_base}={pol3_exp}) but JS={score:.3f}>0.5")

    # Create discriminating position mapping
    disc_pos_map = {pos: (pol2_base, pol3_base) for pos, pol2_base, pol3_base in display_positions}

    # Build the display
    pol2_line = ""
    pol3_line = ""
    consensus_line = ""

    for i in range(len(overall_consensus)):
        if i in disc_pos_map:
            # This is a discriminating position - use the specific bases from discrimination analysis
            pol2_base, pol3_base = disc_pos_map[i]

            pol2_line += pol2_base
            pol3_line += pol3_base
            consensus_line += "|"
        else:
            # Not a discriminating position
            pol2_line += " "
            pol3_line += " "
            consensus_line += overall_consensus[i]

    # Format the compact display
    result = f"Pol II:  {pol2_line}\n"
    result += f"         {consensus_line}\n"
    result += f"Pol III: {pol3_line}"

    # Add warnings if any
    if warnings:
        result += "\n\n" + "\n".join(warnings)

    return result

def analyze_promoter_file(input_file: str, config_file: str = "config.yaml") -> None:
    """Main analysis function"""

    # Load RNA classification config
    rna_config = load_rna_config(config_file)

    # Read sequences
    sequences = []
    seq_ids = []

    try:
        for record in SeqIO.parse(input_file, 'fasta'):
            sequences.append(str(record.seq))
            seq_ids.append(record.id)
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    if not sequences:
        print("No sequences found in input file", file=sys.stderr)
        sys.exit(1)

    # Extract promoter regions
    promoter_seqs, start_pos, end_pos = extract_promoter_regions(sequences, seq_ids)

    # Filter out sequences without promoter regions
    valid_indices = [i for i, seq in enumerate(promoter_seqs) if seq]

    if not valid_indices:
        print("No promoter regions found (no uppercase sequences)", file=sys.stderr)
        sys.exit(1)

    valid_promoters = [promoter_seqs[i] for i in valid_indices]
    valid_ids = [seq_ids[i] for i in valid_indices]

    # Classify sequences
    pol2_promoters = []
    pol3_promoters = []
    pol2_ids = []
    pol3_ids = []

    for i, seq_id in enumerate(valid_ids):
        classification = classify_sequence(seq_id, rna_config)
        if classification == 'pol2':
            pol2_promoters.append(valid_promoters[i])
            pol2_ids.append(seq_id)
        elif classification == 'pol3':
            pol3_promoters.append(valid_promoters[i])
            pol3_ids.append(seq_id)

    # Calculate discriminating positions first to determine consensus strategy
    pol2_freqs = calculate_nucleotide_frequencies(pol2_promoters)
    pol3_freqs = calculate_nucleotide_frequencies(pol3_promoters)
    disc_positions = calculate_position_discrimination(pol2_freqs, pol3_freqs)

    # Create set of discriminatory position indices for quick lookup
    disc_pos_set = {pos for pos, _, _, _ in disc_positions}

    # Generate overall consensus using different strategies per position
    overall_consensus = generate_position_wise_consensus(valid_promoters, disc_pos_set, threshold=0.2)

    # Generate group-specific consensuses
    pol2_consensus = ""
    pol3_consensus = ""
    discriminating_info = "No group-specific analysis available"

    if pol2_promoters and pol3_promoters:
        pol2_consensus = generate_consensus_sequence(pol2_promoters, threshold=0.2)
        pol3_consensus = generate_consensus_sequence(pol3_promoters, threshold=0.2)

        # disc_positions already calculated above
        discriminating_info = format_discriminating_positions(disc_positions)
    else:
        # No discriminating positions if we don't have both groups
        disc_positions = []
        overall_consensus = generate_consensus_sequence(valid_promoters, threshold=0.2)

    # Output results
    print(f"PROMOTER ANALYSIS SUMMARY")
    print(f"========================")
    print(f"Total sequences: {len(valid_promoters)} (Pol II: {len(pol2_promoters)}, Pol III: {len(pol3_promoters)})")
    print(f"Promoter length: {len(overall_consensus)} bp")
    print()

    if pol2_promoters and pol3_promoters:
        # Use compact display for polymerase comparison
        compact_display = create_compact_consensus_display(overall_consensus, pol2_consensus, pol3_consensus, disc_positions)
        print(compact_display)
        print()
        print(f"Key differences: {discriminating_info}")
    else:
        print(f"Overall consensus: {overall_consensus}")
        if pol2_promoters:
            print(f"Pol II only: {pol2_consensus}")
        elif pol3_promoters:
            print(f"Pol III only: {pol3_consensus}")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze promoter consensus sequences and identify Pol II/III differences'
    )

    parser.add_argument('input_file',
                       help='Input FASTA file with aligned promoter sequences (annotated.fa)')
    parser.add_argument('-c', '--config', default='config.yaml',
                       help='Configuration file for RNA classification (default: config.yaml)')

    args = parser.parse_args()

    analyze_promoter_file(args.input_file, args.config)

if __name__ == "__main__":
    main()