#!/usr/bin/env python3
"""
Extract features from GFF file and retrieve sequences from genome with optional extensions.

Usage:
    python3 extract_gff_features.py -g genome.fna.gz -f annotations.gff -t rRNA --extend 50
"""

import sys
import gzip
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO


def load_genome(genome_file: str) -> Dict[str, str]:
    """Load genome sequences into dictionary"""
    sequences = {}

    print(f"Loading genome from {genome_file}...", file=sys.stderr)

    # Handle gzipped or plain files
    if genome_file.endswith('.gz'):
        handle = gzip.open(genome_file, 'rt')
    else:
        handle = open(genome_file, 'r')

    try:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = str(record.seq).upper()
        print(f"Loaded {len(sequences)} sequences", file=sys.stderr)
    finally:
        handle.close()

    return sequences


def parse_gff_for_feature(gff_file: str, feature_type: str, exclude_pseudogenes: bool = True) -> List[Tuple[str, int, int, str, str]]:
    """
    Parse GFF file and extract features of specified type.
    Returns list of (seqid, start, end, strand, attributes)
    """
    features = []
    pseudogene_count = 0

    print(f"Parsing GFF file for {feature_type} features...", file=sys.stderr)

    with open(gff_file, 'r') as f:
        for line in f:
            # Skip comments and empty lines
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid = parts[0]
            feature = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            # Filter by feature type
            if feature == feature_type:
                # Check if this is a pseudogene
                if exclude_pseudogenes:
                    # Check for pseudo=true in attributes
                    if 'pseudo=true' in attributes or 'pseudogene' in attributes.lower():
                        pseudogene_count += 1
                        continue

                features.append((seqid, start, end, strand, attributes))

    print(f"Found {len(features)} {feature_type} features", file=sys.stderr)
    if exclude_pseudogenes and pseudogene_count > 0:
        print(f"Excluded {pseudogene_count} pseudogenes", file=sys.stderr)
    return features


def extract_sequence(sequences: Dict[str, str], seqid: str, start: int, end: int,
                     strand: str, extend: int = 0) -> str:
    """
    Extract sequence from genome with optional extension.
    GFF coordinates are 1-based, inclusive.
    """
    if seqid not in sequences:
        print(f"Warning: Sequence {seqid} not found in genome", file=sys.stderr)
        return ""

    seq = sequences[seqid]
    seq_len = len(seq)

    # Extend coordinates
    ext_start = max(1, start - extend)
    ext_end = min(seq_len, end + extend)

    # Convert to 0-based for Python indexing
    extracted = seq[ext_start - 1:ext_end]

    # Reverse complement if on minus strand
    if strand == '-':
        complement = str.maketrans('ACGT', 'TGCA')
        extracted = extracted.translate(complement)[::-1]

    return extracted


def parse_attribute_value(attributes: str, key: str) -> str:
    """Extract value for a given key from GFF attributes"""
    for attr in attributes.split(';'):
        if '=' in attr:
            k, v = attr.split('=', 1)
            if k == key:
                return v
    return ""


def main():
    parser = argparse.ArgumentParser(
        description='Extract features from GFF and retrieve sequences with extensions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract rRNA with 50nt extensions
  python3 extract_gff_features.py -g genome.fna.gz -f annotations.gff -t rRNA --extend 50

  # Extract tRNA without extensions
  python3 extract_gff_features.py -g genome.fna.gz -f annotations.gff -t tRNA

  # Save to file
  python3 extract_gff_features.py -g genome.fna.gz -f annotations.gff -t rRNA --extend 50 -o output.fa
        """
    )

    parser.add_argument('-g', '--genome', required=True,
                       help='Genome FASTA file (.fna or .fna.gz)')
    parser.add_argument('-f', '--gff', required=True,
                       help='GFF annotation file')
    parser.add_argument('-t', '--type', required=True,
                       help='Feature type to extract (e.g., rRNA, tRNA, CDS)')
    parser.add_argument('--extend', type=int, default=0,
                       help='Number of nucleotides to extend on each side (default: 0)')
    parser.add_argument('--include-pseudogenes', action='store_true',
                       help='Include pseudogenes (default: exclude pseudogenes)')
    parser.add_argument('-o', '--output',
                       help='Output FASTA file (default: stdout)')

    args = parser.parse_args()

    # Check input files exist
    if not Path(args.genome).exists():
        print(f"Error: Genome file not found: {args.genome}", file=sys.stderr)
        sys.exit(1)

    if not Path(args.gff).exists():
        print(f"Error: GFF file not found: {args.gff}", file=sys.stderr)
        sys.exit(1)

    # Load genome
    sequences = load_genome(args.genome)

    # Parse GFF
    exclude_pseudogenes = not args.include_pseudogenes
    features = parse_gff_for_feature(args.gff, args.type, exclude_pseudogenes)

    if not features:
        print(f"No {args.type} features found in GFF file", file=sys.stderr)
        sys.exit(0)

    # Determine output destination
    if args.output:
        try:
            outfile = open(args.output, 'w')
        except IOError as e:
            print(f"Error opening output file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        outfile = sys.stdout

    # Extract sequences
    print(f"Extracting sequences with {args.extend}nt extensions...", file=sys.stderr)
    extracted_count = 0

    try:
        for i, (seqid, start, end, strand, attributes) in enumerate(features, 1):
            # Extract sequence
            seq = extract_sequence(sequences, seqid, start, end, strand, args.extend)

            if seq:
                # Try to get a meaningful ID from attributes
                feature_id = parse_attribute_value(attributes, 'ID')
                gene_name = parse_attribute_value(attributes, 'gene')
                product = parse_attribute_value(attributes, 'product')

                # Build header
                header_parts = [f"{args.type}_{i}"]
                if feature_id:
                    header_parts.append(feature_id)
                if gene_name:
                    header_parts.append(gene_name)

                header = "|".join(header_parts)

                # Add location info
                ext_indicator = f"+{args.extend}nt" if args.extend > 0 else ""
                location = f"{seqid}:{start-args.extend}-{end+args.extend}({strand}){ext_indicator}"

                # Write FASTA (unwrapped)
                outfile.write(f">{header} {location}")
                if product:
                    outfile.write(f" {product}")
                outfile.write("\n")
                outfile.write(f"{seq}\n")
                extracted_count += 1

        print(f"Extracted {extracted_count} sequences", file=sys.stderr)

    finally:
        if args.output:
            outfile.close()


if __name__ == "__main__":
    main()
