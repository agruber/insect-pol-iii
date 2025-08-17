#!/usr/bin/env python3
"""
Genome downloader script - reads filtered TSV data from stdin
Generates download commands for species data received via pipe
"""

import sys
import pandas as pd
import argparse
from datetime import datetime
import io

def clean_species_name(species_name):
    """Convert species name to filesystem-safe format"""
    return species_name.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '_')

def generate_download_commands(df, output_dir="genomes"):
    """Generate download commands for filtered dataframe"""
    commands = []
    
    for _, row in df.iterrows():
        species = clean_species_name(row.iloc[0])  # Species name (first column)
        url = row.iloc[2]  # Download URL (third column)
        
        # Extract filename from URL
        filename = url.split('/')[-1] + "_genomic.fna.gz"
        
        # Create directory structure
        species_dir = f"{output_dir}/{species}"
        logs_dir = f"{species_dir}/logs"
        
        # Download commands
        commands.extend([
            f"mkdir -p {logs_dir}",
            f"wget -O {species_dir}/{filename} {url}/{filename}",
            f"ln -sf {filename} {species_dir}/genome.fna.gz",
            f"echo 'Downloaded from: {url}' > {logs_dir}/download.txt",
            f"echo 'Download date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}' >> {logs_dir}/download.txt",
            f"echo 'Original species name: {row.iloc[0]}' >> {logs_dir}/download.txt"
        ])
    
    return commands

def main():
    parser = argparse.ArgumentParser(description='Generate genome download commands from TSV input')
    parser.add_argument('-o', '--output-dir', default='genomes', help='Output directory for genomes (default: genomes)')
    
    args = parser.parse_args()
    
    # Read TSV data from stdin
    if sys.stdin.isatty():
        print("Error: This script expects TSV data from stdin. Use with filter_species.py", file=sys.stderr)
        print("Example: ./filter_species.py data/genomes_annotated.tsv 'Drosophila melanogaster' | ./download_genomes.py", file=sys.stderr)
        sys.exit(1)
    
    # Read all input
    input_data = sys.stdin.read().strip()
    if not input_data:
        print("Error: No input data received", file=sys.stderr)
        sys.exit(1)
    
    # Parse TSV data
    df = pd.read_csv(io.StringIO(input_data), sep='\t', header=None)
    
    if df.empty:
        print("Error: No valid data received", file=sys.stderr)
        sys.exit(1)
    
    # Generate download commands
    commands = generate_download_commands(df, args.output_dir)
    
    # Print commands
    for cmd in commands:
        print(cmd)

if __name__ == "__main__":
    main()