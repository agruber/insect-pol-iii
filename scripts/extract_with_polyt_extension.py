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
    Extract sequence and extend until poly-T termination signal is found.
    First checks if sequence already has sufficient poly-T (4+ Ts). If not,
    searches downstream/upstream for the nearest poly-T signal and extends to include it.

    Args:
        feature: GFF3 feature with coordinates
        genome_seq: Full genome sequence for the chromosome
        max_extension: Maximum number of bases to extend beyond original end (default: 1000)

    Returns:
        Tuple of (sequence, was_extended, new_start, new_end)
    """
    seq_length = len(genome_seq)

    # Extract original sequence
    start_idx = feature.start - 1  # Convert to 0-based
    end_idx = feature.end  # Already 1-based, so this is correct for slicing

    # Validate original coordinates
    if start_idx < 0 or end_idx > seq_length:
        print(f"Warning: Invalid coordinates {feature.start}-{feature.end} for sequence length {seq_length}", file=sys.stderr)
        return None, False, feature.start, feature.end

    # Extract original sequence
    original_seq = genome_seq[start_idx:end_idx]

    # Apply strand orientation to original sequence
    if feature.strand == '-':
        original_seq = original_seq.reverse_complement()

    # Convert to string and check if it already has sufficient poly-T
    original_str = str(original_seq).upper()

    # Check if sequence already has poly-T termination (4+ consecutive Ts at end)
    if has_polyt_termination(original_str):
        # Already has poly-T, but check if we can extend to get more Ts
        trailing_t_match = re.search(r'T+$', original_str)
        if trailing_t_match:
            if feature.strand == '+':
                # For + strand, check downstream for more Ts
                max_possible = min(max_extension, seq_length - end_idx)
                if max_possible > 0:
                    downstream_seq = str(genome_seq[end_idx:end_idx + max_possible]).upper()
                    leading_t_match = re.match(r'^T+', downstream_seq)
                    if leading_t_match:
                        extension_length = len(leading_t_match.group())
                        new_end = feature.end + extension_length
                        extended_seq = genome_seq[start_idx:end_idx + extension_length]
                        added_sequence = str(genome_seq[end_idx:end_idx + extension_length])
                        print(f"  Extended {feature.seqid}:{feature.start}-{feature.end}({feature.strand}) by {extension_length}bp to {feature.start}-{new_end} (additional Ts found)", file=sys.stderr)
                        print(f"    Added sequence: {added_sequence}", file=sys.stderr)
                        return extended_seq, True, feature.start, new_end
            else:
                # For - strand, check upstream for more As
                max_possible = min(max_extension, start_idx)
                if max_possible > 0:
                    upstream_seq = str(genome_seq[start_idx - max_possible:start_idx]).upper()
                    trailing_a_match = re.search(r'A+$', upstream_seq)
                    if trailing_a_match:
                        extension_length = len(trailing_a_match.group())
                        new_start = feature.start - extension_length
                        extended_region = genome_seq[start_idx - extension_length:end_idx]
                        extended_seq = extended_region.reverse_complement()
                        added_region = genome_seq[start_idx - extension_length:start_idx]
                        added_sequence = str(added_region.reverse_complement())
                        print(f"  Extended {feature.seqid}:{feature.start}-{feature.end}({feature.strand}) by {extension_length}bp to {new_start}-{feature.end} (additional As found upstream)", file=sys.stderr)
                        print(f"    Added sequence: {added_sequence}", file=sys.stderr)
                        return extended_seq, True, new_start, feature.end

        # Already has poly-T and no additional extension possible
        return original_seq, False, feature.start, feature.end

    # Sequence doesn't have poly-T termination, search for it
    print(f"  Searching for poly-T termination for {feature.seqid}:{feature.start}-{feature.end}({feature.strand})", file=sys.stderr)

    if feature.strand == '+':
        # For + strand, search downstream for TTTTT pattern (5+ Ts for extension)
        max_search = min(max_extension, seq_length - end_idx)
        if max_search > 0:
            search_region = str(genome_seq[end_idx:end_idx + max_search]).upper()
            # Look for first occurrence of 5+ consecutive Ts
            polyt_match = re.search(r'T{5,}', search_region)
            if polyt_match:
                # Extend to include the poly-T signal
                extension_to_polyt = polyt_match.end()  # Include the entire poly-T
                new_end = feature.end + extension_to_polyt
                extended_seq = genome_seq[start_idx:end_idx + extension_to_polyt]
                added_sequence = str(genome_seq[end_idx:end_idx + extension_to_polyt])
                print(f"  Extended {feature.seqid}:{feature.start}-{feature.end}({feature.strand}) by {extension_to_polyt}bp to {feature.start}-{new_end} (found 5+ poly-T at +{polyt_match.start()})", file=sys.stderr)
                print(f"    Added sequence: {added_sequence}", file=sys.stderr)
                return extended_seq, True, feature.start, new_end
    else:
        # For - strand, search upstream for AAAAA pattern (becomes TTTTT after reverse complement)
        max_search = min(max_extension, start_idx)
        if max_search > 0:
            search_region = str(genome_seq[start_idx - max_search:start_idx]).upper()
            # Look for last occurrence of 5+ consecutive As (will be first poly-T when reverse complemented)
            polya_matches = list(re.finditer(r'A{5,}', search_region))
            if polya_matches:
                # Take the last match (closest to our sequence)
                polya_match = polya_matches[-1]
                # Extend to include from the start of the poly-A
                extension_to_polya = max_search - polya_match.start()
                new_start = feature.start - extension_to_polya
                extended_region = genome_seq[start_idx - extension_to_polya:end_idx]
                extended_seq = extended_region.reverse_complement()
                added_region = genome_seq[start_idx - extension_to_polya:start_idx]
                added_sequence = str(added_region.reverse_complement())
                print(f"  Extended {feature.seqid}:{feature.start}-{feature.end}({feature.strand}) by {extension_to_polya}bp to {new_start}-{feature.end} (found 5+ poly-A at -{polya_match.start()})", file=sys.stderr)
                print(f"    Added sequence: {added_sequence}", file=sys.stderr)
                return extended_seq, True, new_start, feature.end

    # No 5+ poly-T termination found within search distance
    print(f"  No 5+ poly-T termination found within {max_extension}bp for {feature.seqid}:{feature.start}-{feature.end}({feature.strand})", file=sys.stderr)
    return original_seq, False, feature.start, feature.end

def create_fasta_header(feature, seq_length, species_name, feature_counter, new_start=None, new_end=None):
    """Create simple FASTA header: >Species_name-N|type|seqid|start|end|strand"""
    start = new_start if new_start is not None else feature.start
    end = new_end if new_end is not None else feature.end
    return f"{species_name}-{feature_counter}|{feature.type}|{feature.seqid}|{start}|{end}|{feature.strand}"

def process_gff3_with_extension(gff3_input, genome_dict, output_handle, species_name, incomplete_file=None, extended_file=None, max_extension=1000):
    """Process GFF3 file and extract sequences with poly-T extension"""
    features_processed = 0
    sequences_extracted = 0
    incomplete_sequences = []
    extended_sequences = []

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
        sequence, was_extended, new_start, new_end = extend_sequence_to_polyt(feature, genome_seq, max_extension)

        if sequence is not None:
            # Create FASTA header with potentially extended coordinates
            header = create_fasta_header(feature, len(sequence), species_name, sequences_extracted + 1, new_start, new_end)

            # Track if sequence was extended
            if was_extended:
                seq_id = header.split('|')[0].lstrip('>')
                extended_sequences.append(seq_id)

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

    # Write extended sequences to file if any were found
    if extended_sequences and extended_file:
        with open(extended_file, 'w') as f:
            for seq_id in extended_sequences:
                f.write(f"{seq_id}\n")
        print(f"Wrote {len(extended_sequences)} extended sequence IDs to {extended_file}", file=sys.stderr)

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
    parser.add_argument('--extended',
                       help='Output file for extended sequence IDs (sequences that were extended to find poly-T)')
    parser.add_argument('--max-extension', type=int, default=1000,
                       help='Maximum distance to extend searching for poly-T termination signal (default: 1000)')

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
        process_gff3_with_extension(input_source, genome_dict, output_handle, args.species, args.incomplete, args.extended, args.max_extension)
    finally:
        if args.output:
            output_handle.close()

if __name__ == "__main__":
    main()