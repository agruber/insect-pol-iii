#!/usr/bin/env python3
"""
Filter out sequences containing non-standard nucleotides (not A, C, G, T)

Reads FASTA from stdin and writes filtered FASTA to stdout.
Only keeps sequences containing exclusively A, C, G, T nucleotides.
Statistics are written to stderr.
"""

import sys
import re
from Bio import SeqIO

def has_non_standard_nucleotides(sequence: str) -> bool:
    """Check if sequence contains any nucleotides other than A, C, G, T"""
    # Convert to uppercase and check for non-ACGT characters
    sequence_upper = sequence.upper()
    return bool(re.search(r'[^ACGT]', sequence_upper))

def main():
    total_sequences = 0
    filtered_sequences = 0
    kept_sequences = 0

    for record in SeqIO.parse(sys.stdin, "fasta"):
        total_sequences += 1
        sequence_str = str(record.seq)

        if has_non_standard_nucleotides(sequence_str):
            filtered_sequences += 1
            print(f"Filtered sequence {record.id}: contains non-ACGT nucleotides", file=sys.stderr)
        else:
            kept_sequences += 1
            SeqIO.write(record, sys.stdout, "fasta")

    # Report statistics
    print(f"Non-standard nucleotide filtering statistics:", file=sys.stderr)
    print(f"  Total sequences processed: {total_sequences}", file=sys.stderr)
    print(f"  Sequences with non-ACGT nucleotides: {filtered_sequences}", file=sys.stderr)
    print(f"  Sequences kept: {kept_sequences}", file=sys.stderr)
    if total_sequences > 0:
        print(f"  Retention rate: {kept_sequences/total_sequences:.1%}", file=sys.stderr)

if __name__ == "__main__":
    main()