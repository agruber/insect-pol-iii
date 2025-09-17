#!/usr/bin/env python3
"""
Score sequences using cmsearch with caching.
Scores sequences against data/models/noe.cm and caches results.
"""

import sys
import argparse
import subprocess
import tempfile
import os
import hashlib
from pathlib import Path
from Bio import SeqIO


def get_sequence_hash(sequence):
    """Get MD5 hash of sequence for caching."""
    return hashlib.md5(sequence.encode()).hexdigest()


def load_cache(cache_file):
    """Load existing cache from file."""
    cache = {}
    if os.path.exists(cache_file):
        with open(cache_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and '\t' in line:
                    seq_hash, evalue = line.split('\t', 1)
                    cache[seq_hash] = evalue
    return cache


def save_cache(cache_file, cache):
    """Save cache to file."""
    with open(cache_file, 'w') as f:
        for seq_hash, evalue in cache.items():
            f.write(f"{seq_hash}\t{evalue}\n")


def run_cmsearch(sequence, model_file):
    """Run cmsearch on a single sequence and return e-value."""
    # Create temporary FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_fa:
        tmp_fa.write(f">temp_seq\n{sequence}\n")
        tmp_fa_path = tmp_fa.name

    try:
        # Create temporary tblout file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tblout', delete=False) as tmp_tbl:
            tmp_tbl_path = tmp_tbl.name

        # Run cmsearch with separate tblout file
        cmd = ['cmsearch', '--tblout', tmp_tbl_path, model_file, tmp_fa_path]
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Warning: cmsearch failed with return code {result.returncode}: {result.stderr}", file=sys.stderr)
            return "NA"

        # Parse tblout output for e-value
        with open(tmp_tbl_path, 'r') as tblout_f:
            for line in tblout_f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue

                # tblout format (tab/space-separated): target name, accession, query name, accession, mdl, mdl from, mdl to, seq from, seq to, strand, trunc, pass, gc, bias, score, E-value, inc, description
                fields = line.split()
                if len(fields) >= 16:
                    try:
                        evalue_str = fields[15].strip()
                        # Convert scientific notation to float then back to string to standardize
                        evalue = float(evalue_str)
                        return str(evalue)
                    except (ValueError, IndexError):
                        continue

        # No significant hits found
        return "NA"

    finally:
        # Clean up temporary files
        os.unlink(tmp_fa_path)
        try:
            os.unlink(tmp_tbl_path)
        except:
            pass


def main():
    parser = argparse.ArgumentParser(description='Score sequences using cmsearch with caching')
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('-m', '--model', default='data/models/noe.cm',
                       help='Covariance model file (default: data/models/noe.cm)')
    parser.add_argument('-o', '--output', help='Output scores file (default: same dir as input with .noe.scores suffix)')
    parser.add_argument('-c', '--cache', help='Cache file (default: cache/cmsearch_cache.tsv)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')

    args = parser.parse_args()

    # Set default output file
    if not args.output:
        fasta_path = Path(args.fasta_file)
        args.output = fasta_path.parent / f"{fasta_path.stem}.noe.scores"

    # Set default cache file
    if not args.cache:
        cache_dir = Path("cache")
        cache_dir.mkdir(exist_ok=True)
        args.cache = cache_dir / "cmsearch_cache.tsv"

    # Check if model file exists
    if not os.path.exists(args.model):
        print(f"Error: Model file {args.model} not found", file=sys.stderr)
        sys.exit(1)

    # Check if cmsearch is available
    try:
        subprocess.run(['cmsearch', '-h'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: cmsearch not found. Please install INFERNAL.", file=sys.stderr)
        sys.exit(1)

    # Load cache
    cache = load_cache(args.cache)
    if args.verbose:
        print(f"Loaded {len(cache)} entries from cache", file=sys.stderr)

    # Process sequences
    results = []
    cache_hits = 0
    cache_misses = 0

    for record in SeqIO.parse(args.fasta_file, 'fasta'):
        sequence = str(record.seq).upper()
        seq_hash = get_sequence_hash(sequence)

        if seq_hash in cache:
            evalue = cache[seq_hash]
            cache_hits += 1
            if args.verbose:
                print(f"Cache hit for {record.id}: {evalue}", file=sys.stderr)
        else:
            if args.verbose:
                print(f"Running cmsearch for {record.id}...", file=sys.stderr)
            evalue = run_cmsearch(sequence, args.model)
            cache[seq_hash] = evalue
            cache_misses += 1
            if args.verbose:
                print(f"cmsearch result for {record.id}: {evalue}", file=sys.stderr)

        results.append((record.id, evalue))

    # Save updated cache
    save_cache(args.cache, cache)

    # Write results
    with open(args.output, 'w') as f:
        f.write("Sequence_ID\tE_value\n")
        for seq_id, evalue in results:
            f.write(f"{seq_id}\t{evalue}\n")

    if args.verbose:
        print(f"Results written to {args.output}", file=sys.stderr)
        print(f"Cache statistics: {cache_hits} hits, {cache_misses} misses", file=sys.stderr)
        print(f"Cache saved to {args.cache}", file=sys.stderr)


if __name__ == '__main__':
    main()