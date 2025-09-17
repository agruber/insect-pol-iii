#!/usr/bin/env python3
"""
Find optimal promoter window length by analyzing conservation patterns
around the core k-mer motif.
"""

import argparse
import numpy as np
from Bio import SeqIO
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import yaml
import math
from scipy import stats
import warnings
warnings.filterwarnings('ignore')


def load_rna_config(config_file: str = "config.yaml") -> Dict:
    """Load RNA type configuration"""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config.get('rna_polymerases', {})


def classify_sequence(seq_id: str, rna_config: Dict) -> str:
    """Classify sequence as pol2 or pol3 based on RNA type"""
    # Extract the RNA type from the sequence ID
    # Format: species-number|rna_type|chromosome|start|end|strand
    parts = seq_id.split('|')
    if len(parts) >= 2:
        rna_type = parts[1]

        # Check against pol2 and pol3 types
        pol2_types = rna_config.get('pol2', {}).get('types', [])
        pol3_types = rna_config.get('pol3', {}).get('types', [])

        if rna_type in pol2_types:
            return 'pol2'
        elif rna_type in pol3_types:
            return 'pol3'

    return 'unknown'


def find_kmer_in_sequence(sequence: str, kmer_variants: List[str]) -> Tuple[int, int, str]:
    """Find k-mer position in sequence, return (start, end, matched_kmer)"""
    from itertools import product

    # Expand IUPAC codes in k-mer
    iupac_map = {
        'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
        'M': ['A', 'C'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['C', 'G'],
        'Y': ['C', 'T'], 'K': ['G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'],
        'D': ['A', 'G', 'T'], 'B': ['C', 'G', 'T'], 'N': ['A', 'C', 'G', 'T']
    }

    expanded_kmers = set()
    for kmer in kmer_variants:
        # Expand IUPAC codes
        expanded = ['']
        for base in kmer:
            if base.upper() in iupac_map:
                expanded = [seq + b for seq in expanded for b in iupac_map[base.upper()]]
            else:
                expanded = [seq + base.upper() for seq in expanded]
        expanded_kmers.update(expanded)

    # Search for k-mers in sequence
    seq_upper = sequence.upper()
    for kmer in expanded_kmers:
        pos = seq_upper.find(kmer)
        if pos != -1:
            return pos, pos + len(kmer), kmer

    return -1, -1, ""

def parse_aligned_sequences(fasta_file: str, kmer_variants: List[str], max_extension: int = 18) -> Dict[str, Dict]:
    """
    Parse aligned sequences and find k-mer positions within them.
    Extract extended regions around the k-mer positions.
    """
    sequences = {'all': [], 'pol2': [], 'pol3': []}

    for record in SeqIO.parse(fasta_file, 'fasta'):
        seq = str(record.seq)

        # Find k-mer position in the sequence
        kmer_start, kmer_end, matched_kmer = find_kmer_in_sequence(seq, kmer_variants)

        if kmer_start == -1:
            continue  # Skip sequences without k-mer

        # Extract extended region around k-mer
        left_ext = max(0, kmer_start - max_extension)
        right_ext = min(len(seq), kmer_end + max_extension)

        extended_seq = seq[left_ext:right_ext].upper()

        # Store with metadata
        seq_data = {
            'id': record.id,
            'sequence': extended_seq,
            'core_start': kmer_start - left_ext,
            'core_end': kmer_end - left_ext,
            'matched_kmer': matched_kmer,
            'original_seq': seq,
            'original_kmer_pos': (kmer_start, kmer_end)
        }

        sequences['all'].append(seq_data)

    return sequences


def classify_sequences_by_polymerase(sequences: Dict, rna_config: Dict, verbose: bool = False) -> Dict:
    """Classify sequences by RNA polymerase type"""
    for seq_data in sequences['all']:
        pol_type = classify_sequence(seq_data['id'], rna_config)
        if pol_type == 'pol2':
            sequences['pol2'].append(seq_data)
        elif pol_type == 'pol3':
            sequences['pol3'].append(seq_data)

    return sequences


def calculate_position_frequency_matrix(sequences: List[Dict]) -> np.ndarray:
    """Calculate position frequency matrix"""
    if not sequences:
        return np.array([])

    max_len = max(len(s['sequence']) for s in sequences)
    pfm = np.zeros((4, max_len))
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    for seq_data in sequences:
        seq = seq_data['sequence']
        for i, base in enumerate(seq):
            if base in base_to_idx:
                pfm[base_to_idx[base], i] += 1

    # Normalize
    col_sums = pfm.sum(axis=0)
    col_sums[col_sums == 0] = 1  # Avoid division by zero
    pfm = pfm / col_sums

    return pfm


def pfm_to_consensus(pfm: np.ndarray, left_bound: int = 0, right_bound: int = None) -> str:
    """Convert PFM to consensus sequence"""
    if pfm.size == 0:
        return ""

    if right_bound is None:
        right_bound = pfm.shape[1]

    bases = ['A', 'C', 'G', 'T']
    consensus = ""

    for i in range(left_bound, min(right_bound, pfm.shape[1])):
        max_idx = np.argmax(pfm[:, i])
        max_freq = pfm[max_idx, i]

        # Use IUPAC codes for ambiguous positions
        if max_freq >= 0.75:
            consensus += bases[max_idx]
        elif max_freq >= 0.5:
            consensus += bases[max_idx].lower()
        else:
            # Check for 2-way ambiguity
            sorted_indices = np.argsort(pfm[:, i])[::-1]
            if pfm[sorted_indices[0], i] + pfm[sorted_indices[1], i] >= 0.75:
                b1, b2 = sorted([bases[sorted_indices[0]], bases[sorted_indices[1]]])
                iupac = {'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K'}
                consensus += iupac.get(b1+b2, 'N')
            else:
                consensus += 'N'

    return consensus


def find_discriminatory_positions(pfm_pol2: np.ndarray, pfm_pol3: np.ndarray,
                                 left_bound: int, right_bound: int,
                                 threshold: float = 0.3) -> List[int]:
    """Find positions where Pol II and Pol III differ significantly

    Returns list of positions (relative to left_bound) that are discriminatory
    """
    discriminatory = []

    for i in range(left_bound, min(right_bound, min(pfm_pol2.shape[1], pfm_pol3.shape[1]))):
        # Calculate Jensen-Shannon divergence for this position
        pol2_dist = pfm_pol2[:, i]
        pol3_dist = pfm_pol3[:, i]

        # Skip if either distribution is all zeros
        if np.sum(pol2_dist) == 0 or np.sum(pol3_dist) == 0:
            continue

        # Calculate JS divergence
        m = 0.5 * (pol2_dist + pol3_dist)

        # KL divergence calculations with small epsilon to avoid log(0)
        epsilon = 1e-10
        pol2_dist = pol2_dist + epsilon
        pol3_dist = pol3_dist + epsilon
        m = m + epsilon

        kl_pol2_m = np.sum(pol2_dist * np.log2(pol2_dist / m))
        kl_pol3_m = np.sum(pol3_dist * np.log2(pol3_dist / m))

        js_div = 0.5 * kl_pol2_m + 0.5 * kl_pol3_m

        # If JS divergence exceeds threshold, mark as discriminatory
        # But only if both distributions have some information (not all uniform)
        pol2_entropy = -np.sum(pol2_dist * np.log2(pol2_dist + epsilon))
        pol3_entropy = -np.sum(pol3_dist * np.log2(pol3_dist + epsilon))

        # Only mark as discriminatory if both have reasonable information content
        if js_div >= threshold and pol2_entropy < 1.8 and pol3_entropy < 1.8:
            discriminatory.append(i - left_bound)

    return discriminatory


def enhance_consensus_with_discrimination(consensus: str, discriminatory_positions: List[int],
                                        core_start_rel: int, core_end_rel: int) -> str:
    """Add markers to consensus to show discriminatory positions

    Uses:
    - UPPERCASE: k-mer region (core)
    - [X]: discriminatory position outside k-mer
    - {X}: discriminatory position inside k-mer
    """
    if not consensus:
        return consensus

    enhanced = list(consensus)

    for pos in discriminatory_positions:
        if pos < len(enhanced):
            char = enhanced[pos]

            # Skip marking 'N' positions as discriminatory (they're just uninformative)
            if char.upper() == 'N':
                continue

            # Check if position is in k-mer region
            if core_start_rel <= pos < core_end_rel:
                # Discriminatory position in k-mer: use {X}
                enhanced[pos] = f"{{{char}}}"
            else:
                # Discriminatory position outside k-mer: use [X]
                enhanced[pos] = f"[{char}]"

    return ''.join(enhanced)


def calculate_information_content(pfm: np.ndarray) -> np.ndarray:
    """Calculate information content at each position"""
    ic = np.zeros(pfm.shape[1])

    for i in range(pfm.shape[1]):
        # Shannon entropy
        entropy = 0
        for j in range(4):
            if pfm[j, i] > 0:
                entropy -= pfm[j, i] * np.log2(pfm[j, i])

        # Information content (max 2 bits for DNA)
        ic[i] = 2 - entropy

    return ic


def calculate_relative_entropy(pfm: np.ndarray, background: np.ndarray = None) -> np.ndarray:
    """Calculate relative entropy (KL divergence) against background"""
    if background is None:
        background = np.array([0.25, 0.25, 0.25, 0.25])

    re = np.zeros(pfm.shape[1])

    for i in range(pfm.shape[1]):
        for j in range(4):
            if pfm[j, i] > 0:
                re[i] += pfm[j, i] * np.log2(pfm[j, i] / background[j])

    return re


def randomization_test(sequences: List[Dict], position: int, n_iterations: int = 1000) -> float:
    """Perform randomization test for position significance"""
    if not sequences or position >= len(sequences[0]['sequence']):
        return 1.0

    # Get observed frequencies
    obs_counts = defaultdict(int)
    for seq_data in sequences:
        if position < len(seq_data['sequence']):
            obs_counts[seq_data['sequence'][position]] += 1

    # Calculate observed chi-square
    total = sum(obs_counts.values())
    expected = total / 4
    obs_chi2 = sum((count - expected) ** 2 / expected for count in obs_counts.values())

    # Randomization
    higher_count = 0
    for _ in range(n_iterations):
        rand_counts = defaultdict(int)
        for seq_data in sequences:
            if position < len(seq_data['sequence']):
                # Random base
                rand_counts[np.random.choice(['A', 'C', 'G', 'T'])] += 1

        rand_chi2 = sum((count - expected) ** 2 / expected for count in rand_counts.values())
        if rand_chi2 >= obs_chi2:
            higher_count += 1

    return higher_count / n_iterations


def trim_non_informative_ends(ic_pol2: np.ndarray, ic_pol3: np.ndarray,
                              left_bound: int, right_bound: int,
                              threshold: float = 0.3, verbose: bool = False) -> Tuple[int, int]:
    """Trim non-informative columns from the ends of the window

    Remove positions from ends where both Pol II and Pol III have low information content
    """
    if verbose:
        print(f"\nTrimming non-informative ends from window {left_bound}-{right_bound}")
        print(f"  Trimming threshold: {threshold}")

    # Trim from left end
    new_left = left_bound
    for i in range(left_bound, right_bound):
        pol2_ic = ic_pol2[i] if i < len(ic_pol2) else 0
        pol3_ic = ic_pol3[i] if i < len(ic_pol3) else 0

        # Keep position if either has decent information content
        if pol2_ic >= threshold or pol3_ic >= threshold:
            new_left = i
            break

        if verbose:
            print(f"  Trimming left position {i}: Pol II IC={pol2_ic:.3f}, Pol III IC={pol3_ic:.3f}")

    # Trim from right end
    new_right = right_bound
    for i in range(right_bound - 1, new_left - 1, -1):
        pol2_ic = ic_pol2[i] if i < len(ic_pol2) else 0
        pol3_ic = ic_pol3[i] if i < len(ic_pol3) else 0

        # Keep position if either has decent information content
        # Use higher threshold (0.4) to be more aggressive in trimming
        if pol2_ic >= 0.4 or pol3_ic >= 0.4:
            new_right = i + 1
            break

        if verbose:
            print(f"  Trimming right position {i}: Pol II IC={pol2_ic:.3f}, Pol III IC={pol3_ic:.3f}")

    if verbose:
        print(f"  Final trimmed window: {new_left}-{new_right} (was {left_bound}-{right_bound})")

    return new_left, new_right


def find_unified_conservation_boundary(ic_pol2: np.ndarray, ic_pol3: np.ndarray,
                                      core_start: int, core_end: int,
                                      direction: str, threshold: float = 0.75,
                                      window_size: int = 3, verbose: bool = False) -> int:
    """Find conservation boundary using unified Pol II/III evaluation

    Accept a position if either Pol II OR Pol III exceeds the threshold
    """
    if verbose:
        print(f"Finding {direction} boundary from core region {core_start}-{core_end} (unified Pol II/III)")
        print(f"  Threshold: {threshold}, Window size: {window_size}")

    if direction == 'left':
        # Move left from core
        for i in range(core_start - 1, -1, -1):
            if i >= window_size:
                pol2_window_mean = np.mean(ic_pol2[i-window_size:i]) if len(ic_pol2) > i else 0
                pol3_window_mean = np.mean(ic_pol3[i-window_size:i]) if len(ic_pol3) > i else 0

                # Accept if either Pol II OR Pol III exceeds threshold
                max_conservation = max(pol2_window_mean, pol3_window_mean)

                if verbose:
                    print(f"  Position {i}: Pol II={pol2_window_mean:.3f}, Pol III={pol3_window_mean:.3f}, max={max_conservation:.3f}")

                if max_conservation < threshold:
                    if verbose:
                        print(f"  Left boundary found at position {i + 1}")
                    return i + 1
        if verbose:
            print(f"  Left boundary reached sequence start (0)")
        return 0

    else:  # right
        max_len = min(len(ic_pol2), len(ic_pol3))
        for i in range(core_end, max_len):
            if i + window_size <= max_len:
                pol2_window_mean = np.mean(ic_pol2[i:i+window_size])
                pol3_window_mean = np.mean(ic_pol3[i:i+window_size])

                # Accept if either Pol II OR Pol III exceeds threshold
                max_conservation = max(pol2_window_mean, pol3_window_mean)

                if verbose:
                    print(f"  Position {i}: Pol II={pol2_window_mean:.3f}, Pol III={pol3_window_mean:.3f}, max={max_conservation:.3f}")

                if max_conservation < threshold:
                    if verbose:
                        print(f"  Right boundary found at position {i}")
                    return i
        if verbose:
            print(f"  Right boundary reached sequence end ({max_len})")
        return max_len


def find_conservation_boundary(ic: np.ndarray, core_start: int, core_end: int,
                              direction: str, threshold: float = 0.75,
                              window_size: int = 3, verbose: bool = False) -> int:
    """Find conservation boundary using sliding window"""
    if verbose:
        print(f"Finding {direction} boundary from core region {core_start}-{core_end}")
        print(f"  Threshold: {threshold}, Window size: {window_size}")

    if direction == 'left':
        # Move left from core
        for i in range(core_start - 1, -1, -1):
            if i >= window_size:
                window_mean = np.mean(ic[i-window_size:i])
                if verbose:
                    print(f"  Position {i}: window_mean = {window_mean:.3f}")
                if window_mean < threshold:
                    if verbose:
                        print(f"  Left boundary found at position {i + 1}")
                    return i + 1
        if verbose:
            print(f"  Left boundary reached sequence start (0)")
        return 0

    else:  # right
        # Move right from core
        for i in range(core_end, len(ic)):
            if i + window_size <= len(ic):
                window_mean = np.mean(ic[i:i+window_size])
                if verbose:
                    print(f"  Position {i}: window_mean = {window_mean:.3f}")
                if window_mean < threshold:
                    if verbose:
                        print(f"  Right boundary found at position {i}")
                    return i
        if verbose:
            print(f"  Right boundary reached sequence end ({len(ic)})")
        return len(ic)


def bootstrap_confidence_interval(sequences: List[Dict], n_bootstrap: int = 100) -> Tuple[int, int]:
    """Calculate bootstrap confidence interval for window size"""
    window_sizes = []

    for _ in range(n_bootstrap):
        # Sample with replacement
        sample_indices = np.random.choice(len(sequences), len(sequences), replace=True)
        sample_seqs = [sequences[i] for i in sample_indices]

        # Calculate window for sample
        pfm = calculate_position_frequency_matrix(sample_seqs)
        ic = calculate_information_content(pfm)

        if len(sample_seqs) > 0:
            core_start = sample_seqs[0]['core_start']
            core_end = sample_seqs[0]['core_end']

            left = find_conservation_boundary(ic, core_start, core_end, 'left')
            right = find_conservation_boundary(ic, core_start, core_end, 'right')

            window_sizes.append(right - left)

    if window_sizes:
        return np.percentile(window_sizes, 2.5), np.percentile(window_sizes, 97.5)
    return 0, 0


def test_polymerase_difference(pol2_seqs: List[Dict], pol3_seqs: List[Dict],
                               n_permutations: int = 1000) -> float:
    """Test if Pol II and Pol III have significantly different window sizes"""
    if not pol2_seqs or not pol3_seqs:
        return 1.0

    # Calculate observed windows
    pol2_pfm = calculate_position_frequency_matrix(pol2_seqs)
    pol2_ic = calculate_information_content(pol2_pfm)
    pol3_pfm = calculate_position_frequency_matrix(pol3_seqs)
    pol3_ic = calculate_information_content(pol3_pfm)

    if len(pol2_seqs) > 0 and len(pol3_seqs) > 0:
        pol2_core_start = pol2_seqs[0]['core_start']
        pol2_core_end = pol2_seqs[0]['core_end']
        pol3_core_start = pol3_seqs[0]['core_start']
        pol3_core_end = pol3_seqs[0]['core_end']

        pol2_left = find_conservation_boundary(pol2_ic, pol2_core_start, pol2_core_end, 'left')
        pol2_right = find_conservation_boundary(pol2_ic, pol2_core_start, pol2_core_end, 'right')
        pol3_left = find_conservation_boundary(pol3_ic, pol3_core_start, pol3_core_end, 'left')
        pol3_right = find_conservation_boundary(pol3_ic, pol3_core_start, pol3_core_end, 'right')

        obs_diff = abs((pol2_right - pol2_left) - (pol3_right - pol3_left))
    else:
        return 1.0

    # Permutation test
    all_seqs = pol2_seqs + pol3_seqs
    n_pol2 = len(pol2_seqs)
    greater_count = 0

    for _ in range(n_permutations):
        np.random.shuffle(all_seqs)
        perm_pol2 = all_seqs[:n_pol2]
        perm_pol3 = all_seqs[n_pol2:]

        if perm_pol2 and perm_pol3:
            perm_pol2_pfm = calculate_position_frequency_matrix(perm_pol2)
            perm_pol2_ic = calculate_information_content(perm_pol2_pfm)
            perm_pol3_pfm = calculate_position_frequency_matrix(perm_pol3)
            perm_pol3_ic = calculate_information_content(perm_pol3_pfm)

            p2_left = find_conservation_boundary(perm_pol2_ic, pol2_core_start, pol2_core_end, 'left')
            p2_right = find_conservation_boundary(perm_pol2_ic, pol2_core_start, pol2_core_end, 'right')
            p3_left = find_conservation_boundary(perm_pol3_ic, pol3_core_start, pol3_core_end, 'left')
            p3_right = find_conservation_boundary(perm_pol3_ic, pol3_core_start, pol3_core_end, 'right')

            perm_diff = abs((p2_right - p2_left) - (p3_right - p3_left))
            if perm_diff >= obs_diff:
                greater_count += 1

    return greater_count / n_permutations


def analyze_optimal_window(sequences: Dict, rna_config: Dict, threshold: float = 0.75, verbose: bool = False) -> Dict:
    """Main analysis to find optimal promoter window using unified Pol II/III approach"""
    results = {}

    # First ensure we have both Pol II and Pol III sequences
    if not sequences['pol2'] or not sequences['pol3']:
        if verbose:
            print("Warning: Missing Pol II or Pol III sequences, falling back to overall analysis")
        # Fallback to original approach if missing data
        return analyze_optimal_window_fallback(sequences, rna_config, threshold, verbose)

    if verbose:
        print(f"\nUnified analysis with {len(sequences['pol2'])} Pol II and {len(sequences['pol3'])} Pol III sequences...")

    # Calculate information content for both polymerases
    pfm_pol2 = calculate_position_frequency_matrix(sequences['pol2'])
    ic_pol2 = calculate_information_content(pfm_pol2)
    pfm_pol3 = calculate_position_frequency_matrix(sequences['pol3'])
    ic_pol3 = calculate_information_content(pfm_pol3)

    core_start = sequences['all'][0]['core_start']
    core_end = sequences['all'][0]['core_end']

    if verbose:
        print(f"Core k-mer region: positions {core_start}-{core_end}")
        print(f"Pol II information content range: {ic_pol2.min():.3f} - {ic_pol2.max():.3f}")
        print(f"Pol III information content range: {ic_pol3.min():.3f} - {ic_pol3.max():.3f}")

    # Use unified boundary detection with configurable threshold
    left_bound = find_unified_conservation_boundary(ic_pol2, ic_pol3, core_start, core_end, 'left',
                                                   threshold=threshold, verbose=verbose)
    right_bound = find_unified_conservation_boundary(ic_pol2, ic_pol3, core_start, core_end, 'right',
                                                    threshold=threshold, verbose=verbose)

    # First cap at maximum window length of 25nt
    initial_length = right_bound - left_bound
    if initial_length > 25:
        if verbose:
            print(f"\nInitial window length {initial_length} exceeds maximum of 25nt, capping to 25")
        # Adjust the right boundary to maintain the max length
        right_bound = left_bound + 25
        if verbose:
            print(f"Capped window boundaries: {left_bound} to {right_bound}")

    # Now trim non-informative ends from the capped window
    trimmed_left, trimmed_right = trim_non_informative_ends(ic_pol2, ic_pol3, left_bound, right_bound, verbose=verbose)

    # Adjust to promoter coordinates (negative = upstream)
    promoter_start = trimmed_left - core_end
    promoter_length = trimmed_right - trimmed_left
    # Calculate ps (promoter shift) - nucleotides upstream of k-mer start where promoter begins
    ps_value = core_start - trimmed_left

    if verbose:
        print(f"Final unified window boundaries: {trimmed_left} to {trimmed_right}")
        print(f"Unified promoter window: start={promoter_start}, length={promoter_length}")
        print(f"Unified calculated ps (promoter shift): {ps_value}")

    # Store unified results
    results['unified'] = {
        'start': promoter_start,
        'length': promoter_length,
        'ps_value': ps_value,
        'core_position': f"{core_start-trimmed_left}-{core_end-trimmed_left}"
    }

    # Store individual polymerase results for comparison
    ci_low_pol2, ci_high_pol2 = bootstrap_confidence_interval(sequences['pol2'])
    ci_low_pol3, ci_high_pol3 = bootstrap_confidence_interval(sequences['pol3'])

    results['pol2'] = {
        'start': promoter_start,
        'length': promoter_length,
        'ps_value': ps_value,
        'conservation_profile': ic_pol2[trimmed_left:trimmed_right].tolist() if trimmed_right <= len(ic_pol2) else ic_pol2[trimmed_left:].tolist(),
        'bootstrap_ci': [ci_low_pol2, ci_high_pol2],
        'n_sequences': len(sequences['pol2'])
    }

    results['pol3'] = {
        'start': promoter_start,
        'length': promoter_length,
        'ps_value': ps_value,
        'conservation_profile': ic_pol3[trimmed_left:trimmed_right].tolist() if trimmed_right <= len(ic_pol3) else ic_pol3[trimmed_left:].tolist(),
        'bootstrap_ci': [ci_low_pol3, ci_high_pol3],
        'n_sequences': len(sequences['pol3'])
    }

    # Test for difference between polymerases (using original individual analysis)
    if verbose:
        print(f"\nTesting for differences between Pol II and Pol III...")

    p_value = test_polymerase_difference(sequences['pol2'], sequences['pol3'])
    results['statistical_support'] = {
        'different_windows_pvalue': p_value,
        'significantly_different': p_value < 0.05
    }

    return results


def analyze_optimal_window_fallback(sequences: Dict, rna_config: Dict, threshold: float = 0.75, verbose: bool = False) -> Dict:
    """Fallback to original analysis when Pol II or III data is missing"""
    results = {}

    # Analyze overall if available
    if sequences['all']:
        if verbose:
            print(f"\nFallback: Analyzing overall window from {len(sequences['all'])} sequences...")

        pfm_all = calculate_position_frequency_matrix(sequences['all'])
        ic_all = calculate_information_content(pfm_all)

        core_start = sequences['all'][0]['core_start']
        core_end = sequences['all'][0]['core_end']

        if verbose:
            print(f"Core k-mer region: positions {core_start}-{core_end}")
            print(f"Information content range: {ic_all.min():.3f} - {ic_all.max():.3f}")

        left_bound = find_conservation_boundary(ic_all, core_start, core_end, 'left', verbose=verbose)
        right_bound = find_conservation_boundary(ic_all, core_start, core_end, 'right', verbose=verbose)

        promoter_start = left_bound - core_end
        promoter_length = right_bound - left_bound
        ps_value = core_start - left_bound

        # Cap at maximum window length of 25nt
        if promoter_length > 25:
            if verbose:
                print(f"Window length {promoter_length} exceeds maximum of 25nt, capping to 25")
            right_bound = left_bound + 25
            promoter_length = 25
            ps_value = core_start - left_bound

        if verbose:
            print(f"Window boundaries: {left_bound} to {right_bound}")
            print(f"Promoter window: start={promoter_start}, length={promoter_length}")
            print(f"Calculated ps (promoter shift): {ps_value}")

        results['overall'] = {
            'start': promoter_start,
            'length': promoter_length,
            'ps_value': ps_value,
            'core_position': f"{core_start-left_bound}-{core_end-left_bound}"
        }

    # Individual polymerase analysis for comparison
    for pol_type in ['pol2', 'pol3']:
        if sequences[pol_type]:
            pfm = calculate_position_frequency_matrix(sequences[pol_type])
            ic = calculate_information_content(pfm)

            core_start = sequences[pol_type][0]['core_start']
            core_end = sequences[pol_type][0]['core_end']

            left_bound = find_conservation_boundary(ic, core_start, core_end, 'left')
            right_bound = find_conservation_boundary(ic, core_start, core_end, 'right')

            promoter_start = left_bound - core_end
            promoter_length = right_bound - left_bound
            ps_value = core_start - left_bound

            if promoter_length > 25:
                right_bound = left_bound + 25
                promoter_length = 25
                ps_value = core_start - left_bound

            ci_low, ci_high = bootstrap_confidence_interval(sequences[pol_type])

            results[pol_type] = {
                'start': promoter_start,
                'length': promoter_length,
                'ps_value': ps_value,
                'conservation_profile': ic[left_bound:right_bound].tolist(),
                'bootstrap_ci': [ci_low, ci_high],
                'n_sequences': len(sequences[pol_type])
            }

    if sequences['pol2'] and sequences['pol3']:
        p_value = test_polymerase_difference(sequences['pol2'], sequences['pol3'])
        results['statistical_support'] = {
            'different_windows_pvalue': p_value,
            'significantly_different': p_value < 0.05
        }

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Find optimal promoter window length based on conservation patterns'
    )
    parser.add_argument('input_fasta', help='Aligned sequences from regular annotation')
    parser.add_argument('-c', '--config', default='config.yaml', help='Configuration file')
    parser.add_argument('--kmer', required=True, help='K-mer to search for (supports IUPAC codes)')
    parser.add_argument('-o', '--output', help='Output file for results')
    parser.add_argument('-t', '--threshold', type=float, default=0.75,
                       help='Conservation threshold for boundary detection (default 0.75, range 0-2)')
    parser.add_argument('-e', '--extension', type=int, default=18,
                       help='Maximum extension from k-mer position')
    parser.add_argument('--n-bootstrap', type=int, default=100,
                       help='Number of bootstrap iterations')
    parser.add_argument('--n-permutations', type=int, default=1000,
                       help='Number of permutations for significance testing')
    parser.add_argument('--verbose', action='store_true',
                       help='Show detailed decision process')
    parser.add_argument('--eval', action='store_true',
                       help='Evaluate multiple thresholds and output comparison table')

    args = parser.parse_args()

    # Load configuration
    rna_config = load_rna_config(args.config)

    # If eval mode, run multiple thresholds
    if args.eval:
        print("\nEvaluating conservation thresholds from 0.30 to 1.00 (step 0.05)")
        print("=" * 140)
        print(f"{'Threshold':<12}{'Window':<10}{'PS':<8}{'Start':<10}{'Core':<12}{'Consensus Sequence ([]=[pol2/3 diff], {}={{k-mer diff}})':<50}")
        print("-" * 140)

        for threshold in np.arange(0.30, 1.05, 0.05):
            # Parse sequences
            kmer_variants = [args.kmer]
            sequences = parse_aligned_sequences(args.input_fasta, kmer_variants, args.extension)
            sequences = classify_sequences_by_polymerase(sequences, rna_config)

            # Analyze with this threshold
            results = analyze_optimal_window(sequences, rna_config, threshold=threshold, verbose=False)

            # Calculate consensus for unified sequences and discriminatory positions
            pfm_all = calculate_position_frequency_matrix(sequences['all'])
            pfm_pol2 = calculate_position_frequency_matrix(sequences['pol2'])
            pfm_pol3 = calculate_position_frequency_matrix(sequences['pol3'])

            if 'unified' in results:
                window_length = results['unified']['length']
                ps_value = results['unified']['ps_value']
                start_pos = results['unified']['start']
                core_pos = results['unified']['core_position']

                # Get the actual boundaries used
                core_start = sequences['all'][0]['core_start']
                core_end = sequences['all'][0]['core_end']
                left_bound = core_start - ps_value
                right_bound = left_bound + window_length

                # Generate consensus
                consensus = pfm_to_consensus(pfm_all, left_bound, right_bound)

                # Mark k-mer region in consensus with uppercase
                if len(consensus) >= core_end - left_bound:
                    consensus_list = list(consensus)
                    for i in range(max(0, core_start - left_bound), min(len(consensus), core_end - left_bound)):
                        consensus_list[i] = consensus_list[i].upper()
                    consensus = ''.join(consensus_list)

                # Find discriminatory positions
                if len(sequences['pol2']) > 0 and len(sequences['pol3']) > 0:
                    discriminatory_pos = find_discriminatory_positions(pfm_pol2, pfm_pol3, left_bound, right_bound)
                    core_start_rel = core_start - left_bound
                    core_end_rel = core_end - left_bound
                    consensus = enhance_consensus_with_discrimination(consensus, discriminatory_pos, core_start_rel, core_end_rel)

                # Final cleanup: trim trailing uninformative characters
                consensus = consensus.rstrip('Nn')

                print(f"{threshold:<12.2f}{window_length:<10}{ps_value:<8}{start_pos:<10}{core_pos:<12}{consensus:<50}")
            elif 'overall' in results:
                window_length = results['overall']['length']
                ps_value = results['overall']['ps_value']
                start_pos = results['overall']['start']
                core_pos = results['overall']['core_position']

                # Get the actual boundaries used
                core_start = sequences['all'][0]['core_start']
                core_end = sequences['all'][0]['core_end']
                left_bound = core_start - ps_value
                right_bound = left_bound + window_length

                # Generate consensus
                consensus = pfm_to_consensus(pfm_all, left_bound, right_bound)

                # Mark k-mer region in consensus with uppercase
                if len(consensus) >= core_end - left_bound:
                    consensus_list = list(consensus)
                    for i in range(max(0, core_start - left_bound), min(len(consensus), core_end - left_bound)):
                        consensus_list[i] = consensus_list[i].upper()
                    consensus = ''.join(consensus_list)

                # Find discriminatory positions
                if len(sequences['pol2']) > 0 and len(sequences['pol3']) > 0:
                    discriminatory_pos = find_discriminatory_positions(pfm_pol2, pfm_pol3, left_bound, right_bound)
                    core_start_rel = core_start - left_bound
                    core_end_rel = core_end - left_bound
                    consensus = enhance_consensus_with_discrimination(consensus, discriminatory_pos, core_start_rel, core_end_rel)

                # Final cleanup: trim trailing uninformative characters
                consensus = consensus.rstrip('Nn')

                print(f"{threshold:<12.2f}{window_length:<10}{ps_value:<8}{start_pos:<10}{core_pos:<12}{consensus:<50}")

        print("-" * 140)
        print("\nLegend:")
        print("  UPPERCASE: Highly conserved (>75%) or k-mer region")
        print("  lowercase: Moderately conserved (50-75%)")
        print("  [X]: Pol II/III discriminatory position outside k-mer")
        print("  {X}: Pol II/III discriminatory position inside k-mer")
        print("\nNote: Lower thresholds are more permissive, higher thresholds are stricter")
        print("Recommended range: 0.70-0.80 for balanced results")
        return

    # Parse sequences and find k-mer positions
    kmer_variants = [args.kmer]  # Can be extended to support multiple k-mers
    print(f"Parsing aligned sequences from {args.input_fasta}...")
    print(f"Searching for k-mer: {args.kmer}")
    sequences = parse_aligned_sequences(args.input_fasta, kmer_variants, args.extension)
    print(f"  Found {len(sequences['all'])} sequences with k-mer")

    # Classify by polymerase
    sequences = classify_sequences_by_polymerase(sequences, rna_config, verbose=args.verbose)
    print(f"  Pol II: {len(sequences['pol2'])} sequences")
    print(f"  Pol III: {len(sequences['pol3'])} sequences")

    # Find optimal windows
    print(f"\nAnalyzing optimal promoter windows (threshold={args.threshold})...")
    results = analyze_optimal_window(sequences, rna_config, threshold=args.threshold, verbose=args.verbose)

    # Output results
    output_text = "OPTIMAL PROMOTER WINDOWS\n"
    output_text += "=" * 50 + "\n\n"

    if 'unified' in results:
        output_text += "Unified (Pol II + Pol III):\n"
        output_text += f"  Start position: {results['unified']['start']}\n"
        output_text += f"  Window length: {results['unified']['length']} nt\n"
        output_text += f"  Promoter shift (ps): {results['unified']['ps_value']}\n"
        output_text += f"  Core position: {results['unified']['core_position']}\n\n"
    elif 'overall' in results:
        output_text += "Overall:\n"
        output_text += f"  Start position: {results['overall']['start']}\n"
        output_text += f"  Window length: {results['overall']['length']} nt\n"
        output_text += f"  Promoter shift (ps): {results['overall']['ps_value']}\n"
        output_text += f"  Core position: {results['overall']['core_position']}\n\n"

    if 'pol2' in results:
        output_text += "Pol II (for comparison):\n"
        output_text += f"  Start position: {results['pol2']['start']}\n"
        output_text += f"  Window length: {results['pol2']['length']} nt\n"
        output_text += f"  Promoter shift (ps): {results['pol2']['ps_value']}\n"
        output_text += f"  Bootstrap CI: {results['pol2']['bootstrap_ci']}\n"
        output_text += f"  N sequences: {results['pol2']['n_sequences']}\n\n"

    if 'pol3' in results:
        output_text += "Pol III (for comparison):\n"
        output_text += f"  Start position: {results['pol3']['start']}\n"
        output_text += f"  Window length: {results['pol3']['length']} nt\n"
        output_text += f"  Promoter shift (ps): {results['pol3']['ps_value']}\n"
        output_text += f"  Bootstrap CI: {results['pol3']['bootstrap_ci']}\n"
        output_text += f"  N sequences: {results['pol3']['n_sequences']}\n\n"

    if 'statistical_support' in results:
        output_text += "Statistical Support:\n"
        output_text += f"  Different windows p-value: {results['statistical_support']['different_windows_pvalue']:.4f}\n"
        output_text += f"  Significantly different: {results['statistical_support']['significantly_different']}\n"

    print(output_text)

    if args.output:
        # Save detailed results as YAML
        with open(args.output, 'w') as f:
            yaml.dump(results, f, default_flow_style=False)
        print(f"\nDetailed results saved to {args.output}")


if __name__ == '__main__':
    main()