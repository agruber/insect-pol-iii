#!/usr/bin/env python3
"""
Fast genome screening for Pol II and Pol III promoters using optimized implementations.
Automatically uses the fastest available implementation (Rust > C > Python).
"""

import sys
import gzip
import json
import argparse
import subprocess
from pathlib import Path
from Bio import SeqIO

# Import fallback Python implementation
sys.path.insert(0, str(Path(__file__).parent))
from match_algorithm import MATCHAlgorithm, load_pwm

def detect_best_implementation():
    """Detect the fastest available MATCH implementation"""
    rust_binary = Path("scripts/target/release/match_algorithm_rust")
    c_binary = Path("scripts/match_algorithm_c")

    if rust_binary.exists():
        return "rust", str(rust_binary)
    elif c_binary.exists():
        return "c", str(c_binary)
    else:
        return "python", None

def parse_gff3_attributes(attr_string: str) -> dict:
    """Parse GFF3 attributes string"""
    attributes = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key] = value
    return attributes

def screen_with_external_binary(binary_path: str, impl_name: str,
                               genome_file: Path, pwm_files: dict,
                               cutoff_strategy: str, min_score: float) -> list:
    """Screen using external binary (Rust or C)"""
    predictions = []

    for pwm_type, pwm_file in pwm_files.items():
        if not pwm_file or not Path(pwm_file).exists():
            continue

        print(f"  Running {impl_name} implementation for {pwm_type}...", file=sys.stderr)

        if impl_name == "rust":
            cmd = [binary_path, "-p", pwm_file, "-g", str(genome_file),
                   "-c", cutoff_strategy, "-s", str(min_score)]
        else:  # C implementation
            cmd = [binary_path, "-p", pwm_file, "-g", str(genome_file),
                   "-c", cutoff_strategy, "-s", str(min_score)]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Parse GFF3 output
            for line in result.stdout.split('\n'):
                if line.strip() and not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 9:
                        predictions.append({
                            'seqid': parts[0],
                            'source': f'MATCH_{impl_name.capitalize()}_{pwm_type}',
                            'type': 'promoter',
                            'start': int(parts[3]),
                            'end': int(parts[4]),
                            'score': float(parts[5]),
                            'strand': parts[6],
                            'attributes': parse_gff3_attributes(parts[8])
                        })

        except subprocess.CalledProcessError as e:
            print(f"Error running {impl_name} implementation for {pwm_type}: {e.stderr}", file=sys.stderr)
            # Fall back to Python for this PWM
            python_predictions = screen_with_python(genome_file, {pwm_type: pwm_file},
                                                   cutoff_strategy, min_score)
            predictions.extend(python_predictions)

        except Exception as e:
            print(f"Error with {impl_name} implementation for {pwm_type}: {e}", file=sys.stderr)
            # Fall back to Python for this PWM
            python_predictions = screen_with_python(genome_file, {pwm_type: pwm_file},
                                                   cutoff_strategy, min_score)
            predictions.extend(python_predictions)

    return predictions

def screen_with_python(genome_file: Path, pwm_files: dict,
                      cutoff_strategy: str, min_score: float) -> list:
    """Screen using Python implementation (fallback)"""
    predictions = []

    # Load PWMs
    pwms = {}
    for pwm_type, pwm_file in pwm_files.items():
        if pwm_file and Path(pwm_file).exists():
            try:
                pwms[pwm_type] = load_pwm(pwm_file)
                print(f"  Running Python implementation for {pwm_type}...", file=sys.stderr)
            except Exception as e:
                print(f"Error loading {pwm_type} PWM: {e}", file=sys.stderr)
                continue

    if not pwms:
        return predictions

    # Initialize MATCH algorithms
    matchers = {}
    for pwm_type, pwm in pwms.items():
        matchers[pwm_type] = MATCHAlgorithm(pwm, cutoff_strategy)

    # Open genome file
    if str(genome_file).endswith('.gz'):
        handle = gzip.open(genome_file, 'rt')
    else:
        handle = open(genome_file, 'r')

    # Process each sequence
    for record in SeqIO.parse(handle, 'fasta'):
        seq_id = record.id
        sequence = str(record.seq)
        seq_len = len(sequence)

        for pwm_type, matcher in matchers.items():
            matches = matcher.scan_sequence(sequence)

            for match in matches:
                if match['score'] >= min_score:
                    predictions.append({
                        'seqid': seq_id,
                        'source': f'MATCH_Python_{pwm_type}',
                        'type': 'promoter',
                        'start': match['start'] + 1,  # Convert to 1-based
                        'end': match['end'],
                        'score': match['score'],
                        'strand': match['strand'],
                        'attributes': {
                            'ID': f"{pwm_type}_promoter_{seq_id}_{match['start']}",
                            'polymerase': pwm_type,
                            'css': match['css'],
                            'mss': match['mss'],
                            'sequence': match['sequence']
                        }
                    })

    handle.close()
    return predictions

def screen_genome_fast(genome_file: Path, pol2_pwm_file: str = None, pol3_pwm_file: str = None,
                      cutoff_strategy: str = 'minSum', min_score: float = 0.8) -> list:
    """
    Screen genome using the fastest available implementation.
    """
    # Detect best implementation
    impl_type, binary_path = detect_best_implementation()

    print(f"Using {impl_type.upper()} implementation for genome screening", file=sys.stderr)

    # Prepare PWM files
    pwm_files = {}
    if pol2_pwm_file:
        pwm_files['pol2'] = pol2_pwm_file
    if pol3_pwm_file:
        pwm_files['pol3'] = pol3_pwm_file

    if not pwm_files:
        print("Warning: No PWM files provided", file=sys.stderr)
        return []

    # Screen with the best implementation
    if impl_type in ['rust', 'c']:
        return screen_with_external_binary(binary_path, impl_type, genome_file,
                                         pwm_files, cutoff_strategy, min_score)
    else:
        return screen_with_python(genome_file, pwm_files, cutoff_strategy, min_score)

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

def write_summary(predictions: list, output_file: str = None):
    """Write summary statistics"""
    output = sys.stderr if output_file is None else open(output_file, 'w')

    # Count by implementation and polymerase
    impl_counts = {}
    pol_counts = {'pol2': 0, 'pol3': 0}

    for pred in predictions:
        source = pred['source']
        impl_counts[source] = impl_counts.get(source, 0) + 1

        # Extract polymerase type
        if 'pol2' in source:
            pol_counts['pol2'] += 1
        elif 'pol3' in source:
            pol_counts['pol3'] += 1

    # Write summary
    output.write("=== Fast Promoter Screening Summary ===\n")
    output.write(f"Total predictions: {len(predictions)}\n")
    output.write(f"  Pol II promoters: {pol_counts['pol2']}\n")
    output.write(f"  Pol III promoters: {pol_counts['pol3']}\n")

    if impl_counts:
        output.write(f"\nImplementation breakdown:\n")
        for impl, count in sorted(impl_counts.items()):
            output.write(f"  {impl}: {count}\n")

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
        description='Fast genome screening for Pol II and Pol III promoters'
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
    parser.add_argument('--force-python', action='store_true',
                       help='Force use of Python implementation')

    args = parser.parse_args()

    # Check that at least one PWM is provided
    if not args.pol2_pwm and not args.pol3_pwm:
        print("Error: At least one PWM file must be provided", file=sys.stderr)
        sys.exit(1)

    genome_path = Path(args.genome)
    if not genome_path.exists():
        print(f"Error: Genome file {genome_path} not found", file=sys.stderr)
        sys.exit(1)

    # Screen genome
    if args.force_python:
        print("Forcing Python implementation", file=sys.stderr)
        pwm_files = {}
        if args.pol2_pwm:
            pwm_files['pol2'] = args.pol2_pwm
        if args.pol3_pwm:
            pwm_files['pol3'] = args.pol3_pwm
        predictions = screen_with_python(genome_path, pwm_files, args.cutoff, args.min_score)
    else:
        predictions = screen_genome_fast(genome_path, args.pol2_pwm, args.pol3_pwm,
                                        args.cutoff, args.min_score)

    # Write output (only GFF3 for now since it's most useful)
    write_gff3(predictions, args.output)

    # Write summary if requested
    if args.summary:
        write_summary(predictions, args.summary)
    else:
        write_summary(predictions)

if __name__ == "__main__":
    main()