#!/usr/bin/env python3
"""
Check if there's sufficient ncRNA diversity for annotation.
Requires at least 3 different ncRNA types for both Pol II and Pol III.

Reads FASTA from stdin or file and checks sequence headers for diversity.
Returns exit code 0 if sufficient diversity, 1 if not.
"""

import sys
import argparse
import yaml
from Bio import SeqIO
from pathlib import Path
from collections import defaultdict

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

def check_diversity(fasta_file, config_file: str = "config.yaml", verbose: bool = False) -> tuple:
    """
    Check if there's sufficient ncRNA diversity.
    Returns (has_sufficient_diversity, pol2_types, pol3_types, message)
    """
    polymerase_config = load_polymerase_config(config_file)

    if not polymerase_config:
        return False, set(), set(), "Could not load polymerase configuration"

    pol2_types = set(polymerase_config.get('pol2', {}).get('types', []))
    pol3_types = set(polymerase_config.get('pol3', {}).get('types', []))

    # Count ncRNA types found
    found_pol2 = set()
    found_pol3 = set()
    total_sequences = 0

    # Read sequences
    if fasta_file == '-':
        input_handle = sys.stdin
    else:
        input_handle = open(fasta_file, 'r')

    for record in SeqIO.parse(input_handle, "fasta"):
        total_sequences += 1
        ncrna_type = parse_ncrna_from_header(record.id)

        if ncrna_type:
            if ncrna_type in pol2_types:
                found_pol2.add(ncrna_type)
            elif ncrna_type in pol3_types:
                found_pol3.add(ncrna_type)
            else:
                if verbose:
                    print(f"Warning: Unknown ncRNA type '{ncrna_type}' in {record.id}", file=sys.stderr)

    if fasta_file != '-':
        input_handle.close()

    # Check diversity requirements
    pol2_count = len(found_pol2)
    pol3_count = len(found_pol3)

    has_sufficient = pol2_count >= 3 and pol3_count >= 3

    # Build message
    message_parts = []
    message_parts.append(f"Total sequences: {total_sequences}")
    message_parts.append(f"Pol II ncRNAs: {pol2_count} types ({', '.join(sorted(found_pol2)) if found_pol2 else 'none'})")
    message_parts.append(f"Pol III ncRNAs: {pol3_count} types ({', '.join(sorted(found_pol3)) if found_pol3 else 'none'})")

    if not has_sufficient:
        if pol2_count < 3:
            message_parts.append(f"INSUFFICIENT: Need at least 3 Pol II types, found {pol2_count}")
        if pol3_count < 3:
            message_parts.append(f"INSUFFICIENT: Need at least 3 Pol III types, found {pol3_count}")
    else:
        message_parts.append("SUFFICIENT: Meets diversity requirements")

    message = "\n".join(message_parts)

    return has_sufficient, found_pol2, found_pol3, message

def main():
    parser = argparse.ArgumentParser(
        description='Check ncRNA diversity for annotation requirements'
    )
    parser.add_argument('input', nargs='?', default='-',
                       help='Input FASTA file (default: stdin)')
    parser.add_argument('-c', '--config', default='config.yaml',
                       help='Configuration file (default: config.yaml)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output')
    parser.add_argument('--write-flag',
                       help='Write flag file to this path if insufficient diversity')

    args = parser.parse_args()

    # Check diversity
    has_sufficient, pol2_types, pol3_types, message = check_diversity(
        args.input, args.config, args.verbose
    )

    # Print results
    print(message, file=sys.stderr)

    # Write flag file if requested and insufficient
    if args.write_flag and not has_sufficient:
        flag_path = Path(args.write_flag)
        with open(flag_path, 'w') as f:
            f.write("Insufficient ncRNA diversity for annotation\n")
            f.write(message + "\n")
        print(f"Flag file written to: {flag_path}", file=sys.stderr)

    # Exit with appropriate code
    if has_sufficient:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()