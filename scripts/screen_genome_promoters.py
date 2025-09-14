#!/usr/bin/env python3
"""
Screen genomes for Pol II and Pol III promoters using MATCH algorithm.
Outputs predictions in GFF3 format.
"""

import sys
import gzip
import json
import argparse
from pathlib import Path
from Bio import SeqIO
sys.path.insert(0, str(Path(__file__).parent))
from match_algorithm import MATCHAlgorithm, load_pwm

def screen_genome(genome_file: Path, pol2_pwm: dict, pol3_pwm: dict,
                 cutoff_strategy: str = 'minSum', min_score: float = 0.8) -> list:
    """
    Screen a genome for Pol II and Pol III promoters.

    Args:
        genome_file: Path to genome FASTA file (can be gzipped)
        pol2_pwm: Pol II PWM dictionary
        pol3_pwm: Pol III PWM dictionary
        cutoff_strategy: MATCH cutoff strategy
        min_score: Minimum combined score to report

    Returns:
        List of promoter predictions
    """
    predictions = []

    # Initialize MATCH algorithms
    pol2_matcher = MATCHAlgorithm(pol2_pwm, cutoff_strategy) if pol2_pwm else None
    pol3_matcher = MATCHAlgorithm(pol3_pwm, cutoff_strategy) if pol3_pwm else None

    # Open genome file
    if str(genome_file).endswith('.gz'):
        handle = gzip.open(genome_file, 'rt')
    else:
        handle = open(genome_file, 'r')

    print(f"Screening genome {genome_file}...", file=sys.stderr)

    # Process each sequence
    for record in SeqIO.parse(handle, 'fasta'):
        seq_id = record.id
        sequence = str(record.seq)
        seq_len = len(sequence)

        print(f"  Processing {seq_id} ({seq_len:,} bp)...", file=sys.stderr)

        # Screen with Pol II PWM
        if pol2_matcher:
            pol2_matches = pol2_matcher.scan_sequence(sequence)

            for match in pol2_matches:
                if match['score'] >= min_score:
                    predictions.append({
                        'seqid': seq_id,
                        'source': 'MATCH_pol2',
                        'type': 'promoter',
                        'start': match['start'] + 1,  # Convert to 1-based
                        'end': match['end'],
                        'score': match['score'],
                        'strand': match['strand'],
                        'attributes': {
                            'ID': f"pol2_promoter_{seq_id}_{match['start']}",
                            'polymerase': 'pol2',
                            'css': match['css'],
                            'mss': match['mss'],
                            'sequence': match['sequence']
                        }
                    })

        # Screen with Pol III PWM
        if pol3_matcher:
            pol3_matches = pol3_matcher.scan_sequence(sequence)

            for match in pol3_matches:
                if match['score'] >= min_score:
                    predictions.append({
                        'seqid': seq_id,
                        'source': 'MATCH_pol3',
                        'type': 'promoter',
                        'start': match['start'] + 1,  # Convert to 1-based
                        'end': match['end'],
                        'score': match['score'],
                        'strand': match['strand'],
                        'attributes': {
                            'ID': f"pol3_promoter_{seq_id}_{match['start']}",
                            'polymerase': 'pol3',
                            'css': match['css'],
                            'mss': match['mss'],
                            'sequence': match['sequence']
                        }
                    })

    handle.close()

    print(f"Found {len(predictions)} promoter predictions", file=sys.stderr)

    # Sort by sequence and position
    predictions.sort(key=lambda x: (x['seqid'], x['start']))

    return predictions

def write_gff3(predictions: list, output_file: str = None):
    """Write predictions in GFF3 format"""
    output = sys.stdout if output_file is None else open(output_file, 'w')

    # Write header
    output.write("##gff-version 3\n")

    # Write predictions
    for pred in predictions:
        # Build attributes string
        attrs = []
        for key, value in pred['attributes'].items():
            if key != 'sequence':  # Don't include sequence in GFF
                attrs.append(f"{key}={value}")

        gff_line = "\t".join([
            pred['seqid'],
            pred['source'],
            pred['type'],
            str(pred['start']),
            str(pred['end']),
            f"{pred['score']:.3f}",
            pred['strand'],
            '.',  # phase
            ';'.join(attrs)
        ])

        output.write(gff_line + "\n")

    if output_file:
        output.close()

def write_bed(predictions: list, output_file: str = None):
    """Write predictions in BED format"""
    output = sys.stdout if output_file is None else open(output_file, 'w')

    # Write header
    output.write("track name=promoters description=\"Predicted Pol II and Pol III promoters\"\n")

    # Write predictions
    for pred in predictions:
        # BED is 0-based, half-open
        bed_line = "\t".join([
            pred['seqid'],
            str(pred['start'] - 1),  # Convert to 0-based
            str(pred['end']),
            pred['attributes']['ID'],
            str(int(pred['score'] * 1000)),  # Scale score to 0-1000
            pred['strand']
        ])

        output.write(bed_line + "\n")

    if output_file:
        output.close()

def write_summary(predictions: list, output_file: str = None):
    """Write summary statistics"""
    output = sys.stderr if output_file is None else open(output_file, 'w')

    # Count by polymerase
    pol2_count = sum(1 for p in predictions if p['attributes']['polymerase'] == 'pol2')
    pol3_count = sum(1 for p in predictions if p['attributes']['polymerase'] == 'pol3')

    # Count by sequence
    seq_counts = {}
    for pred in predictions:
        seq_id = pred['seqid']
        pol = pred['attributes']['polymerase']
        if seq_id not in seq_counts:
            seq_counts[seq_id] = {'pol2': 0, 'pol3': 0}
        seq_counts[seq_id][pol] += 1

    # Write summary
    output.write("=== Promoter Screening Summary ===\n")
    output.write(f"Total predictions: {len(predictions)}\n")
    output.write(f"  Pol II promoters: {pol2_count}\n")
    output.write(f"  Pol III promoters: {pol3_count}\n")
    output.write(f"\nSequences with predictions: {len(seq_counts)}\n")

    # Score distribution
    if predictions:
        scores = [p['score'] for p in predictions]
        output.write(f"\nScore statistics:\n")
        output.write(f"  Min: {min(scores):.3f}\n")
        output.write(f"  Max: {max(scores):.3f}\n")
        output.write(f"  Mean: {sum(scores)/len(scores):.3f}\n")

    if output_file:
        output.close()

def main():
    parser = argparse.ArgumentParser(
        description='Screen genomes for Pol II and Pol III promoters'
    )
    parser.add_argument('genome', help='Genome FASTA file (can be gzipped)')
    parser.add_argument('--pol2-pwm', help='Pol II PWM JSON file')
    parser.add_argument('--pol3-pwm', help='Pol III PWM JSON file')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    parser.add_argument('-f', '--format', choices=['gff3', 'bed'],
                       default='gff3', help='Output format (default: gff3)')
    parser.add_argument('-c', '--cutoff', choices=['minFN', 'minFP', 'minSum'],
                       default='minSum', help='Cutoff strategy (default: minSum)')
    parser.add_argument('-s', '--min-score', type=float, default=0.8,
                       help='Minimum score threshold (default: 0.8)')
    parser.add_argument('--summary', help='Write summary to file')

    args = parser.parse_args()

    # Check that at least one PWM is provided
    if not args.pol2_pwm and not args.pol3_pwm:
        print("Error: At least one PWM file must be provided", file=sys.stderr)
        sys.exit(1)

    # Load PWMs
    pol2_pwm = None
    pol3_pwm = None

    if args.pol2_pwm:
        print(f"Loading Pol II PWM from {args.pol2_pwm}...", file=sys.stderr)
        pol2_pwm = load_pwm(args.pol2_pwm)

    if args.pol3_pwm:
        print(f"Loading Pol III PWM from {args.pol3_pwm}...", file=sys.stderr)
        pol3_pwm = load_pwm(args.pol3_pwm)

    # Screen genome
    genome_path = Path(args.genome)
    if not genome_path.exists():
        print(f"Error: Genome file {genome_path} not found", file=sys.stderr)
        sys.exit(1)

    predictions = screen_genome(
        genome_path,
        pol2_pwm,
        pol3_pwm,
        args.cutoff,
        args.min_score
    )

    # Write output
    if args.format == 'gff3':
        write_gff3(predictions, args.output)
    else:
        write_bed(predictions, args.output)

    # Write summary if requested
    if args.summary:
        write_summary(predictions, args.summary)
    else:
        write_summary(predictions)

if __name__ == "__main__":
    main()