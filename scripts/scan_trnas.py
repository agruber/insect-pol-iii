#!/usr/bin/env python3
"""
Generate tRNAscan-SE commands for downloaded genomes
Reads TSV data from stdin and generates bash commands to run tRNAscan-SE
"""

import sys
import os
import gzip
import tempfile
import subprocess
import pandas as pd
from pathlib import Path

def clean_species_name(name):
    """Clean species name by removing problematic characters"""
    return name.replace(' ', '_').replace('(', '').replace(')', '').replace('[', '').replace(']', '')

def run_trnascan(genome_path, output_dir, species_name):
    """
    Run tRNAscan-SE on a genome file
    
    Args:
        genome_path: Path to the genome file (may be gzipped)
        output_dir: Directory to store results
        species_name: Clean species name for logging
    
    Returns:
        bool: True if successful, False otherwise
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Output files
    gff_output = os.path.join(output_dir, "trnascan_results.gff")
    gff_output_gz = gff_output + ".gz"
    log_file = os.path.join(output_dir, "logs", "trnascan.log")
    
    # Create logs directory
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    
    # Check if output already exists
    if os.path.exists(gff_output_gz):
        print(f"tRNAscan results already exist for {species_name}, skipping...")
        return True
    
    try:
        # Handle gzipped genome files
        if genome_path.endswith('.gz'):
            # Create temporary uncompressed file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as temp_file:
                temp_file_path = temp_file.name
                
                print(f"Decompressing {genome_path}...")
                with gzip.open(genome_path, 'rt') as gz_file:
                    temp_file.write(gz_file.read())
        else:
            temp_file_path = genome_path
        
        print(f"Running tRNAscan-SE on {species_name}...")
        
        # Run tRNAscan-SE: -E (eukaryotic), -L (legacy mode), --gff (GFF3 output), -q (quiet)
        with open(log_file, 'w') as log:
            result = subprocess.run([
                "tRNAscan-SE", "-q", "-E", "-L", "--gff", gff_output, temp_file_path
            ], stdout=log, stderr=subprocess.STDOUT, check=True)
        
        # Clean up temporary file if we created one
        if genome_path.endswith('.gz') and os.path.exists(temp_file_path):
            os.unlink(temp_file_path)
        
        # Compress the GFF output
        if os.path.exists(gff_output):
            print(f"Compressing tRNAscan results for {species_name}...")
            with open(gff_output, 'rb') as f_in:
                with gzip.open(gff_output_gz, 'wb') as f_out:
                    f_out.write(f_in.read())
            
            # Remove uncompressed file
            os.unlink(gff_output)
            
            print(f"✓ tRNAscan-SE completed for {species_name}")
            return True
        else:
            print(f"✗ No output generated for {species_name}")
            return False
            
    except subprocess.CalledProcessError as e:
        print(f"✗ tRNAscan-SE failed for {species_name}: {e}")
        return False
    except Exception as e:
        print(f"✗ Error processing {species_name}: {e}")
        return False

def generate_trnascan_commands(df, output_dir="genomes"):
    """
    Generate tRNAscan-SE commands for species in dataframe
    
    Args:
        df: DataFrame with species data
        output_dir: Base directory for genomes
    """
    
    commands = []
    
    for _, row in df.iterrows():
        species = clean_species_name(row.iloc[0])
        species_dir = f"{output_dir}/{species}"
        
        # Look for genome file (try different possible names)
        genome_candidates = [
            f"{species_dir}/genome.fna.gz",
            f"{species_dir}/genome.fna",
        ]
        
        # Also check for actual downloaded genome files
        if os.path.exists(species_dir):
            for file in os.listdir(species_dir):
                if file.endswith(('_genomic.fna.gz', '_genomic.fna')):
                    genome_candidates.append(f"{species_dir}/{file}")
        
        genome_file = None
        for candidate in genome_candidates:
            if os.path.exists(candidate):
                genome_file = candidate
                break
        
        if genome_file:
            # Create command
            log_dir = f"{species_dir}/logs"
            commands.extend([
                f"mkdir -p {log_dir}",
                f"echo 'Running tRNAscan-SE for {species}...' | tee -a {log_dir}/trnascan.log",
                f"python3 -c \"",
                f"import sys; sys.path.append('scripts')",
                f"from scan_trnas import run_trnascan",
                f"success = run_trnascan('{genome_file}', '{species_dir}', '{species}')",
                f"sys.exit(0 if success else 1)",
                f"\" 2>&1 | tee -a {log_dir}/trnascan.log",
                f"echo 'tRNAscan-SE completed for {species}' >> {log_dir}/trnascan.log",
                ""
            ])
        else:
            print(f"Warning: No genome file found for {species}", file=sys.stderr)
            commands.append(f"echo 'Warning: No genome file found for {species}' >&2")
    
    return commands

def main():
    """Main function to process stdin and generate commands"""
    
    if sys.stdin.isatty():
        print("Usage: cat species.tsv | python3 scan_trnas.py", file=sys.stderr)
        print("   or: python3 scan_trnas.py < species.tsv", file=sys.stderr)
        print("\nThis script reads TSV data from stdin and generates tRNAscan-SE commands", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Read TSV data from stdin
        df = pd.read_csv(sys.stdin, sep='\t', header=None)
        
        if df.empty:
            print("No data received from stdin", file=sys.stderr)
            sys.exit(1)
        
        # Generate tRNAscan-SE commands
        commands = generate_trnascan_commands(df)
        
        # Output commands
        for command in commands:
            print(command)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()