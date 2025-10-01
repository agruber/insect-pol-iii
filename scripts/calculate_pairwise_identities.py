#!/usr/bin/env python3
"""
Calculate pairwise sequence identities using EMBOSS needle alignment
Reads FASTA from stdin and outputs TSV to stdout
"""

import sys
import subprocess
import tempfile
import os
from Bio import SeqIO
from itertools import combinations
import re

def run_needle(seq1_file, seq2_file):
    """
    Run EMBOSS needle to align two sequences and extract identity percentage
    """
    # Create temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.needle', delete=False) as tmp_out:
        tmp_out_path = tmp_out.name

    try:
        # Run needle with default parameters
        # -gapopen 10.0 -gapextend 0.5 are defaults
        cmd = [
            'needle',
            '-asequence', seq1_file,
            '-bsequence', seq2_file,
            '-outfile', tmp_out_path,
            '-auto'  # Run in non-interactive mode
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Error running needle: {result.stderr}", file=sys.stderr)
            return None

        # Parse the output to extract identity percentage
        with open(tmp_out_path, 'r') as f:
            content = f.read()

        # Look for the Identity line in needle output
        # Format: # Identity:     123/456 (27.0%) or # Identity:     123/456 ( 8.9%)
        identity_match = re.search(r'# Identity:\s+\d+/\d+\s+\(\s*(\d+\.?\d*)\s*%\)', content)

        if identity_match:
            return float(identity_match.group(1))
        else:
            print(f"Warning: Could not parse identity from needle output", file=sys.stderr)
            return None

    finally:
        # Clean up temporary file
        if os.path.exists(tmp_out_path):
            os.unlink(tmp_out_path)

def calculate_pairwise_identities():
    """
    Main function to calculate all pairwise identities
    """
    # Read sequences from stdin
    sequences = []
    for record in SeqIO.parse(sys.stdin, "fasta"):
        sequences.append(record)

    if len(sequences) < 2:
        print("Error: Need at least 2 sequences for pairwise comparison", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(sequences)} sequences...", file=sys.stderr)

    # Print header
    print("Sequence1\tSequence2\tIdentity")

    # Calculate all pairwise comparisons
    for seq1, seq2 in combinations(sequences, 2):
        # Create temporary files for the two sequences
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp1:
            SeqIO.write(seq1, tmp1, "fasta")
            tmp1_path = tmp1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp2:
            SeqIO.write(seq2, tmp2, "fasta")
            tmp2_path = tmp2.name

        try:
            # Run needle and get identity
            identity = run_needle(tmp1_path, tmp2_path)

            if identity is not None:
                # Extract clean sequence IDs (everything before first |)
                seq1_id = seq1.id.split('|')[0]
                seq2_id = seq2.id.split('|')[0]

                print(f"{seq1_id}\t{seq2_id}\t{identity:.1f}")

        finally:
            # Clean up temporary files
            if os.path.exists(tmp1_path):
                os.unlink(tmp1_path)
            if os.path.exists(tmp2_path):
                os.unlink(tmp2_path)

    print(f"Completed {len(list(combinations(range(len(sequences)), 2)))} pairwise comparisons", file=sys.stderr)

if __name__ == "__main__":
    calculate_pairwise_identities()