#!/usr/bin/env python3
"""
Benchmark script to compare MATCH algorithm implementations:
- Python implementation
- C implementation
- Rust implementation
"""

import time
import subprocess
import os
import sys
from pathlib import Path
import tempfile

def run_command_with_timing(command, description):
    """Run command and measure execution time"""
    print(f"\n=== {description} ===", file=sys.stderr)
    print(f"Command: {command}", file=sys.stderr)

    start_time = time.time()
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        end_time = time.time()

        execution_time = end_time - start_time

        if result.returncode == 0:
            # Count output lines (matches found)
            output_lines = len([line for line in result.stdout.split('\n')
                              if line.strip() and not line.startswith('#')])
            print(f"✓ SUCCESS: {execution_time:.2f} seconds, {output_lines} matches", file=sys.stderr)
            return execution_time, output_lines, True
        else:
            print(f"✗ FAILED: {result.stderr}", file=sys.stderr)
            return execution_time, 0, False

    except Exception as e:
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"✗ ERROR: {e}", file=sys.stderr)
        return execution_time, 0, False

def create_test_data(test_size="small"):
    """Create test PWM and genome files"""

    # Use existing PWM
    pwm_file = "genomes/Drosophila_melanogaster/test_pwm_pol3.json"
    if not os.path.exists(pwm_file):
        print("Error: Test PWM file not found. Run PWM building first.", file=sys.stderr)
        sys.exit(1)

    # Create test genome of appropriate size
    if test_size == "small":
        # 100KB test - extract from existing genome
        genome_file = "/tmp/test_genome_small.fa"
        subprocess.run(f"zcat genomes/Drosophila_melanogaster/genome.fna.gz | head -n 2000 > {genome_file}",
                      shell=True, check=True)
    elif test_size == "medium":
        # 1MB test
        genome_file = "/tmp/test_genome_medium.fa"
        subprocess.run(f"zcat genomes/Drosophila_melanogaster/genome.fna.gz | head -n 20000 > {genome_file}",
                      shell=True, check=True)
    elif test_size == "large":
        # 10MB test
        genome_file = "/tmp/test_genome_large.fa"
        subprocess.run(f"zcat genomes/Drosophila_melanogaster/genome.fna.gz | head -n 200000 > {genome_file}",
                      shell=True, check=True)
    else:
        # Use full genome
        genome_file = "genomes/Drosophila_melanogaster/genome.fna.gz"

    return pwm_file, genome_file

def benchmark_implementations(test_size="small"):
    """Benchmark all three implementations"""

    print(f"Creating test data ({test_size})...", file=sys.stderr)
    pwm_file, genome_file = create_test_data(test_size)

    # Get file size for reference
    file_size = os.path.getsize(genome_file) / (1024 * 1024)  # MB
    print(f"Test genome size: {file_size:.1f} MB", file=sys.stderr)

    results = {}

    # 1. Python implementation
    python_cmd = f"python3 scripts/screen_genome_promoters.py {genome_file} --pol3-pwm {pwm_file} -c minFN -s 0.6"
    time_py, matches_py, success_py = run_command_with_timing(python_cmd, "Python Implementation")
    results['Python'] = {'time': time_py, 'matches': matches_py, 'success': success_py}

    # 2. C implementation
    if os.path.exists("scripts/match_algorithm_c"):
        c_cmd = f"scripts/match_algorithm_c -p {pwm_file} -g {genome_file} -c minFN -s 0.6"
        time_c, matches_c, success_c = run_command_with_timing(c_cmd, "C Implementation")
        results['C'] = {'time': time_c, 'matches': matches_c, 'success': success_c}
    else:
        print("C implementation not found (compile with gcc)", file=sys.stderr)
        results['C'] = {'time': 0, 'matches': 0, 'success': False}

    # 3. Rust implementation
    if os.path.exists("scripts/target/release/match_algorithm_rust"):
        rust_cmd = f"scripts/target/release/match_algorithm_rust -p {pwm_file} -g {genome_file} -c minFN -s 0.6"
        time_rust, matches_rust, success_rust = run_command_with_timing(rust_cmd, "Rust Implementation")
        results['Rust'] = {'time': time_rust, 'matches': matches_rust, 'success': success_rust}
    else:
        print("Rust implementation not found (build with cargo)", file=sys.stderr)
        results['Rust'] = {'time': 0, 'matches': 0, 'success': False}

    # Clean up temporary files
    if genome_file.startswith("/tmp/"):
        os.remove(genome_file)

    return results, file_size

def print_results(results, file_size):
    """Print benchmark results in a nice table"""

    print(f"\n{'='*60}")
    print(f"BENCHMARK RESULTS ({file_size:.1f} MB genome)")
    print(f"{'='*60}")

    # Table header
    print(f"{'Implementation':<15} {'Time (s)':<10} {'Matches':<10} {'Speedup':<10} {'Status':<10}")
    print(f"{'-'*60}")

    # Find baseline (Python) time for speedup calculation
    python_time = results.get('Python', {}).get('time', 0)

    # Sort by time (fastest first)
    sorted_results = sorted(results.items(), key=lambda x: x[1]['time'] if x[1]['success'] else float('inf'))

    for impl, data in sorted_results:
        if data['success']:
            speedup = f"{python_time/data['time']:.1f}x" if data['time'] > 0 and python_time > 0 else "-"
            status = "✓"
        else:
            speedup = "-"
            status = "✗"

        print(f"{impl:<15} {data['time']:<10.2f} {data['matches']:<10} {speedup:<10} {status:<10}")

    # Additional statistics
    print(f"\n{'='*60}")
    print("PERFORMANCE SUMMARY:")

    successful_results = {k: v for k, v in results.items() if v['success']}
    if successful_results:
        fastest = min(successful_results.items(), key=lambda x: x[1]['time'])
        slowest = max(successful_results.items(), key=lambda x: x[1]['time'])

        print(f"  Fastest: {fastest[0]} ({fastest[1]['time']:.2f}s)")
        print(f"  Slowest: {slowest[0]} ({slowest[1]['time']:.2f}s)")

        if fastest[1]['time'] > 0 and slowest[1]['time'] > 0:
            overall_speedup = slowest[1]['time'] / fastest[1]['time']
            print(f"  Overall speedup: {overall_speedup:.1f}x")

        # Throughput (MB/s)
        print(f"\nTHROUGHPUT (MB/s):")
        for impl, data in successful_results.items():
            if data['time'] > 0:
                throughput = file_size / data['time']
                print(f"  {impl}: {throughput:.1f} MB/s")

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark MATCH implementations')
    parser.add_argument('-s', '--size', choices=['small', 'medium', 'large', 'full'],
                       default='small', help='Test size (default: small)')
    parser.add_argument('--run-all', action='store_true',
                       help='Run all sizes for comprehensive benchmark')

    args = parser.parse_args()

    if args.run_all:
        sizes = ['small', 'medium', 'large']
        for size in sizes:
            print(f"\n{'#'*80}", file=sys.stderr)
            print(f"RUNNING {size.upper()} BENCHMARK", file=sys.stderr)
            print(f"{'#'*80}", file=sys.stderr)

            results, file_size = benchmark_implementations(size)
            print_results(results, file_size)
    else:
        results, file_size = benchmark_implementations(args.size)
        print_results(results, file_size)

if __name__ == "__main__":
    main()