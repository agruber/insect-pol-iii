#!/usr/bin/env python3
"""
Extract FASTA sequences from genome file based on GFF3 coordinates with poly-T extension.
Extends sequences until a poly-T termination signal (4+ consecutive Ts) is found.
"""

import sys
import re
from bioseq_lib import (
    load_fasta_sequences, parse_gff3_file,
    write_fasta_sequence, progress_reporter, smart_open,
    create_standard_parser
)

def has_polyt_termination(sequence, min_ts=4):
    """Check if sequence ends with poly-T termination signal"""
    # Convert to uppercase for comparison
    seq_upper = str(sequence).upper()

    # Check if it ends with 4 or more consecutive Ts
    match = re.search(r'T{4,}$', seq_upper)
    return match is not None

def extend_sequence_to_polyt(feature, genome_seq, max_extension=1000):
    """
    Extract sequence and extend until poly-T termination is found.

    Args:
        feature: GFF3 feature with coordinates
        genome_seq: Full genome sequence for the chromosome
        max_extension: Maximum number of bases to extend beyond original end

    Returns:
        Extended sequence or original sequence if poly-T not found
    """
    seq_length = len(genome_seq)

    # Extract original sequence
    start_idx = feature.start - 1  # Convert to 0-based
    end_idx = feature.end  # Already 1-based, so this is correct for slicing

    # Validate original coordinates
    if start_idx < 0 or end_idx > seq_length:
        print(f"Warning: Invalid coordinates {feature.start}-{feature.end} for sequence length {seq_length}", file=sys.stderr)
        return None

    # Extract original sequence
    original_seq = genome_seq[start_idx:end_idx]

    # Apply strand orientation to original sequence
    if feature.strand == '-':
        original_seq = original_seq.reverse_complement()

    # Check if original sequence already has poly-T termination
    if has_polyt_termination(original_seq):
        return original_seq

    # For extension, we need to work in the direction of transcription
    # For + strand: extend towards higher coordinates (3' direction)
    # For - strand: extend towards lower coordinates (5' direction, but 3' of transcript)

    extended_seq = original_seq
    extension_length = 0

    while extension_length < max_extension:
        extension_length += 1

        if feature.strand == '+':
            # Extend towards higher coordinates
            new_end = end_idx + extension_length
            if new_end > seq_length:
                break
            extended_seq = genome_seq[start_idx:new_end]
        else:
            # Extend towards lower coordinates, then reverse complement
            new_start = start_idx - extension_length
            if new_start < 0:
                break
            extended_region = genome_seq[new_start:end_idx]
            extended_seq = extended_region.reverse_complement()

        # Check if we found poly-T termination
        if has_polyt_termination(extended_seq):
            print(f"  Extended {feature.seqid}:{feature.start}-{feature.end}({feature.strand}) by {extension_length}bp to find poly-T termination", file=sys.stderr)
            return extended_seq

    # If we couldn't find poly-T within max_extension, return original sequence
    print(f"  Warning: Could not find poly-T termination within {max_extension}bp for {feature.seqid}:{feature.start}-{feature.end}({feature.strand})", file=sys.stderr)
    return original_seq

def create_fasta_header(feature, seq_length, species_name, feature_counter):
    """Create simple FASTA header: >Species_name-N|type|seqid|start|end|strand"""
    return f"{species_name}-{feature_counter}|{feature.type}|{feature.seqid}|{feature.start}|{feature.end}|{feature.strand}"

def process_gff3_with_extension(gff3_input, genome_dict, output_handle, species_name, incomplete_file=None):
    """Process GFF3 file and extract sequences with poly-T extension"""
    features_processed = 0
    sequences_extracted = 0
    incomplete_sequences = []

    # Handle both filename and file handle
    if hasattr(gff3_input, 'read'):
        # It's a file handle (like stdin)
        gff3_lines = gff3_input.readlines()
        # Create a temporary file-like object from the lines
        import io
        gff3_content = ''.join(gff3_lines)
        temp_file = io.StringIO(gff3_content)

        # Parse features directly from content
        features = []
        for line in gff3_lines:
            line = line.strip()
            if line and not line.startswith('#'):
                from bioseq_lib import GFF3Feature
                try:
                    feature = GFF3Feature.from_line(line)
                    if feature:
                        features.append(feature)
                except:
                    continue
    else:
        # It's a filename
        features = list(parse_gff3_file(gff3_input))

    for feature in features:
        features_processed += 1

        # Check if we have the chromosome/scaffold in genome
        if feature.seqid not in genome_dict:
            print(f"Warning: Sequence '{feature.seqid}' not found in genome", file=sys.stderr)
            continue

        genome_seq = genome_dict[feature.seqid].seq

        # Extract sequence with poly-T extension
        sequence = extend_sequence_to_polyt(feature, genome_seq)

        if sequence is not None:
            # Create FASTA header
            header = create_fasta_header(feature, len(sequence), species_name, sequences_extracted + 1)

            # Check if sequence is incomplete (doesn't end with poly-T)
            if not has_polyt_termination(sequence):
                incomplete_sequences.append(header.split('|')[0].lstrip('>'))  # Just the ID part
                print(f"  Marking as incomplete: {header.split('|')[0].lstrip('>')}", file=sys.stderr)

            # Write sequence
            write_fasta_sequence(output_handle, header, sequence, lowercase=True)
            sequences_extracted += 1

        progress_reporter(features_processed)

    print(f"Total: processed {features_processed} features, extracted {sequences_extracted} sequences",
          file=sys.stderr)

    # Write incomplete sequences to file if any were found
    if incomplete_sequences and incomplete_file:
        with open(incomplete_file, 'w') as f:
            for seq_id in incomplete_sequences:
                f.write(f"{seq_id}\n")
        print(f"Wrote {len(incomplete_sequences)} incomplete sequence IDs to {incomplete_file}", file=sys.stderr)

def main():
    parser = create_standard_parser(
        description='Extract FASTA sequences from genome file with poly-T extension',
        epilog="""
Examples:
  python3 extract_with_polyt_extension.py -g genome.fna.gz -s species_name < features.gff3 > sequences.fasta
        """
    )

    # Add required arguments
    parser.add_argument('-g', '--genome', required=True,
                       help='Input genome FASTA file (supports .gz compression)')
    parser.add_argument('-s', '--species', required=True,
                       help='Species name for FASTA headers (e.g., Drosophila_melanogaster)')
    parser.add_argument('-f', '--features',
                       help='Input GFF3 file (use stdin if not provided)')
    parser.add_argument('-o', '--output',
                       help='Output FASTA file (use stdout if not provided)')
    parser.add_argument('--incomplete',
                       help='Output file for incomplete sequence IDs (sequences without poly-T termination)')

    args = parser.parse_args()

    print("Loading genome sequences...", file=sys.stderr)
    genome_dict = load_fasta_sequences(args.genome)
    print(f"Loaded {len(genome_dict)} sequences from {args.genome}", file=sys.stderr)

    # Handle output
    output_handle = sys.stdout if not args.output else open(args.output, 'w')

    try:
        # Determine input source
        input_source = args.features if args.features else sys.stdin
        # Process GFF3 and extract sequences
        process_gff3_with_extension(input_source, genome_dict, output_handle, args.species, args.incomplete)
    finally:
        if args.output:
            output_handle.close()

if __name__ == "__main__":
    main()