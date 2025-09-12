#!/usr/bin/env python3
"""
Sort FASTA sequences by RNA polymerase type based on config.yaml.
Reads from STDIN, writes to STDOUT.
Pol III sequences come first, then Pol II sequences.
"""

import sys
import yaml
import argparse
from pathlib import Path
from Bio import SeqIO


def load_config(config_path):
    """Load RNA polymerase configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def get_rna_type(record_id):
    """Extract RNA type from sequence ID.
    
    Expected format: Species-N|RNA_TYPE|...
    Example: Drosophila_melanogaster-1|U1|...
    """
    parts = record_id.split('|')
    if len(parts) >= 2:
        return parts[1]
    return None


def sort_sequences(sequences, config, anchor_sort=False):
    """Sort sequences by polymerase type (pol3 first, then pol2).
    
    Within each polymerase group, maintain the order defined in config.yaml.
    Sequences with unrecognized RNA types are placed at the end.
    If anchor_sort is True, non-anchored sequences (no uppercase) go to the end.
    """
    pol2_types = config['rna_polymerases']['pol2']['types']
    pol3_types = config['rna_polymerases']['pol3']['types']
    
    # Create priority map for sorting
    priority = {}
    
    # Pol3 comes first (lower priority numbers)
    for i, rna_type in enumerate(pol3_types):
        priority[rna_type] = i
    
    # Pol2 comes second
    offset = len(pol3_types)
    for i, rna_type in enumerate(pol2_types):
        priority[rna_type] = offset + i
    
    # Sort function
    def sort_key(record):
        # Check if sequence has uppercase letters (is anchored)
        has_uppercase = any(c.isupper() for c in str(record.seq)) if anchor_sort else True
        
        rna_type = get_rna_type(record.id)
        if not has_uppercase:
            # Non-anchored sequences go to the very end
            return (2, 0, record.id)
        elif rna_type and rna_type in priority:
            return (0, priority[rna_type], record.id)  # Known anchored types first
        else:
            return (1, 0, record.id)  # Unknown anchored types second
    
    return sorted(sequences, key=sort_key)


def main():
    parser = argparse.ArgumentParser(
        description='Sort FASTA sequences by RNA polymerase type (pol3 first, then pol2)',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-c', '--config', 
                        default='config.yaml',
                        help='Path to config.yaml file (default: config.yaml)')
    parser.add_argument('-a', '--anchor',
                        action='store_true',
                        help='Only include sequences that have uppercase letters (anchored sequences)')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Print sorting information to stderr')
    
    args = parser.parse_args()
    
    # Check if config file exists
    config_path = Path(args.config)
    if not config_path.exists():
        # Try in parent directory if not found in current
        config_path = Path(__file__).parent.parent / 'config.yaml'
        if not config_path.exists():
            sys.stderr.write(f"Error: Config file not found at {args.config} or {config_path}\n")
            sys.exit(1)
    
    # Load config
    try:
        config = load_config(config_path)
    except Exception as e:
        sys.stderr.write(f"Error loading config file: {e}\n")
        sys.exit(1)
    
    # Read sequences from stdin
    sequences = list(SeqIO.parse(sys.stdin, 'fasta'))
    
    if not sequences:
        sys.stderr.write("Error: No sequences provided\n")
        sys.exit(1)
    
    if args.verbose:
        sys.stderr.write(f"Loaded {len(sequences)} sequences\n")
    
    # Sort sequences (anchor flag affects sort order, not filtering)
    sorted_sequences = sort_sequences(sequences, config, anchor_sort=args.anchor)
    
    if args.verbose:
        # Report sorting
        pol2_types = config['rna_polymerases']['pol2']['types']
        pol3_types = config['rna_polymerases']['pol3']['types']
        
        pol2_count = 0
        pol3_count = 0
        unknown_count = 0
        anchored_count = 0
        non_anchored_count = 0
        
        for seq in sorted_sequences:
            rna_type = get_rna_type(seq.id)
            has_uppercase = any(c.isupper() for c in str(seq.seq))
            
            if has_uppercase:
                anchored_count += 1
            else:
                non_anchored_count += 1
                
            if rna_type in pol3_types:
                pol3_count += 1
            elif rna_type in pol2_types:
                pol2_count += 1
            else:
                unknown_count += 1
        
        sys.stderr.write(f"\nSequences by polymerase:\n")
        sys.stderr.write(f"  Pol III: {pol3_count} sequences\n")
        sys.stderr.write(f"  Pol II:  {pol2_count} sequences\n")
        if unknown_count > 0:
            sys.stderr.write(f"  Unknown: {unknown_count} sequences\n")
        
        if args.anchor:
            sys.stderr.write(f"\nSequences by anchor status:\n")
            sys.stderr.write(f"  Anchored:     {anchored_count} sequences\n")
            sys.stderr.write(f"  Non-anchored: {non_anchored_count} sequences\n")
        
        sys.stderr.write(f"\nOutput order:\n")
        for i, seq in enumerate(sorted_sequences[:10], 1):
            rna_type = get_rna_type(seq.id)
            pol = "Pol III" if rna_type in pol3_types else "Pol II" if rna_type in pol2_types else "Unknown"
            sys.stderr.write(f"  {i:2}. {seq.id[:50]:<50} [{pol}]\n")
        if len(sorted_sequences) > 10:
            sys.stderr.write(f"  ... and {len(sorted_sequences) - 10} more\n")
    
    # Write sorted sequences to stdout without line wrapping
    for record in sorted_sequences:
        sys.stdout.write(f">{record.description}\n")
        sys.stdout.write(f"{str(record.seq)}\n")


if __name__ == '__main__':
    main()