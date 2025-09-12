#!/usr/bin/env python3
"""
Pipeline runner for RNA analysis using configuration from config.yaml
Executes various analysis pipelines based on command and species
"""

import sys
import os
import argparse
import yaml
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

def load_config(config_file: str = "config.yaml") -> Dict:
    """Load configuration from YAML file"""
    if not os.path.exists(config_file):
        print(f"Error: Configuration file '{config_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        print(f"Error loading config file: {e}", file=sys.stderr)
        sys.exit(1)

def get_taxonomy_info(species_name: str) -> Optional[str]:
    """Get taxonomic group for a species (e.g., Insecta for Drosophila)"""
    # For now, we'll assume all species are Insecta
    # This could be enhanced to read from a taxonomy file
    return "Insecta"

def build_evalue_filters(config: Dict) -> List[str]:
    """Build e-value filter arguments from config"""
    evalues = config.get('evalues', {})
    filters = []
    
    # Build gene-specific e-value filters
    for gene, evalue in evalues.items():
        filters.extend(['-g', f'{gene}:{evalue}'])
    
    return filters

def get_excluded_genes(config: Dict, species_name: str) -> List[str]:
    """Get list of genes to exclude based on species taxonomy"""
    exclude_genes = config.get('exclude_genes', {})
    taxonomy = get_taxonomy_info(species_name)
    
    excluded = []
    if taxonomy and taxonomy in exclude_genes:
        excluded = exclude_genes[taxonomy]
    
    return excluded

def check_required_files(species_dir: Path, required_files: List[str]) -> bool:
    """Check if all required files exist"""
    missing = []
    for file_pattern in required_files:
        file_path = species_dir / file_pattern
        # Handle wildcard patterns
        if '*' in file_pattern:
            matching_files = list(species_dir.glob(file_pattern))
            if not matching_files:
                missing.append(file_pattern)
        elif not file_path.exists():
            missing.append(file_pattern)
    
    if missing:
        print(f"Error: Missing required files in {species_dir}:", file=sys.stderr)
        for f in missing:
            print(f"  - {f}", file=sys.stderr)
        return False
    
    return True

def run_upstream_pipeline(species_name: str, config: Dict, upstream_distance: Optional[int] = None, 
                         default_evalue: float = 1e-7, output_file: Optional[str] = None):
    """Run the upstream extraction pipeline for a species"""
    
    species_dir = Path("genomes") / species_name
    
    # Check required files
    required_files = ["*tblout.gz", "genome.tsv", "genome.fna.gz"]
    if not check_required_files(species_dir, required_files):
        sys.exit(1)
    
    # Use upstream distance from config if not specified
    if upstream_distance is None:
        upstream_distance = config.get('upstream_length', 100)
    
    # Build output path
    if output_file is None:
        output_file = species_dir / "upstream.fa"
    
    # Get configuration
    excluded_genes = get_excluded_genes(config, species_name)
    evalue_filters = build_evalue_filters(config)
    
    # Build pipeline command
    pipeline_parts = []
    
    # 1. Extract from tblout files
    tblout_files = list(species_dir.glob("*tblout.gz"))
    if len(tblout_files) == 1:
        pipeline_parts.append(f"zcat {tblout_files[0]}")
    else:
        # If multiple files, concatenate them
        pipeline_parts.append(f"zcat {species_dir}/*tblout.gz")
    
    # 2. Convert to GFF3
    pipeline_parts.append("./scripts/tblout_to_gff3.py")
    
    # 3. Filter out excluded genes
    if excluded_genes:
        filter_args = ' '.join([f'-r {gene}' for gene in excluded_genes])
        pipeline_parts.append(f"./scripts/filter_genes_gff3.py {filter_args}")
    
    # 4. Apply e-value filters
    filter_cmd = f"./scripts/filter_gff3.py -e {default_evalue}"
    if evalue_filters:
        filter_cmd += ' ' + ' '.join(evalue_filters)
    pipeline_parts.append(filter_cmd)
    
    # 5. Extract upstream regions
    pipeline_parts.append(f"./scripts/extract_upstream_regions.py -u {upstream_distance} -l {species_dir}/genome.tsv")
    
    # 6. Convert to FASTA
    pipeline_parts.append(f"./scripts/gff3_to_fasta.py -g {species_dir}/genome.fna.gz -s {species_name}")
    
    # 7. Sort sequences
    pipeline_parts.append("./scripts/sort_seqs.py")
    
    # Join pipeline with pipes
    full_command = " | ".join(pipeline_parts)
    
    # Add output redirection
    full_command += f" > {output_file}"
    
    # Print the command for transparency
    print(f"Executing pipeline for {species_name}:", file=sys.stderr)
    print(f"  Command: {full_command}", file=sys.stderr)
    if excluded_genes:
        print(f"  Excluding genes: {', '.join(excluded_genes)}", file=sys.stderr)
    print(f"  Output: {output_file}", file=sys.stderr)
    
    # Execute the pipeline
    try:
        result = subprocess.run(full_command, shell=True, capture_output=False, text=True)
        if result.returncode != 0:
            print(f"Error: Pipeline failed with exit code {result.returncode}", file=sys.stderr)
            sys.exit(1)
        print(f"Pipeline completed successfully. Output written to {output_file}", file=sys.stderr)
    except Exception as e:
        print(f"Error executing pipeline: {e}", file=sys.stderr)
        sys.exit(1)

def run_filter_pipeline(species_name: str, config: Dict, input_file: Optional[str] = None,
                       output_file: Optional[str] = None, default_evalue: float = 1e-7):
    """Run a simple filtering pipeline on existing tblout or GFF3 files"""
    
    species_dir = Path("genomes") / species_name
    
    # Determine input
    if input_file is None:
        # Look for tblout files
        tblout_files = list(species_dir.glob("*tblout.gz"))
        if not tblout_files:
            print(f"Error: No tblout.gz files found in {species_dir}", file=sys.stderr)
            sys.exit(1)
        input_cmd = f"zcat {' '.join(str(f) for f in tblout_files)} | ./scripts/tblout_to_gff3.py"
    else:
        # Use provided input file
        if input_file.endswith('.gz'):
            input_cmd = f"zcat {input_file}"
        else:
            input_cmd = f"cat {input_file}"
    
    # Build output path
    if output_file is None:
        output_file = species_dir / "filtered.gff3"
    
    # Get configuration
    excluded_genes = get_excluded_genes(config, species_name)
    evalue_filters = build_evalue_filters(config)
    
    # Build pipeline
    pipeline_parts = [input_cmd]
    
    # Filter out excluded genes
    if excluded_genes:
        filter_args = ' '.join([f'-r {gene}' for gene in excluded_genes])
        pipeline_parts.append(f"./scripts/filter_genes_gff3.py {filter_args}")
    
    # Apply e-value filters
    filter_cmd = f"./scripts/filter_gff3.py -e {default_evalue}"
    if evalue_filters:
        filter_cmd += ' ' + ' '.join(evalue_filters)
    pipeline_parts.append(filter_cmd)
    
    # Join pipeline
    full_command = " | ".join(pipeline_parts) + f" > {output_file}"
    
    print(f"Executing filter pipeline for {species_name}:", file=sys.stderr)
    print(f"  Command: {full_command}", file=sys.stderr)
    
    # Execute
    try:
        result = subprocess.run(full_command, shell=True, capture_output=False, text=True)
        if result.returncode != 0:
            print(f"Error: Pipeline failed with exit code {result.returncode}", file=sys.stderr)
            sys.exit(1)
        print(f"Filter pipeline completed. Output: {output_file}", file=sys.stderr)
    except Exception as e:
        print(f"Error executing pipeline: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Run RNA analysis pipelines using configuration from config.yaml',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Commands:
  upstream    Extract upstream regions and generate FASTA file
  filter      Apply filtering to tblout or GFF3 files

Examples:
  # Extract upstream regions for Drosophila melanogaster
  python3 run_pipeline.py upstream Drosophila_melanogaster
  
  # Extract with custom upstream distance
  python3 run_pipeline.py upstream Drosophila_melanogaster -u 200
  
  # Extract with custom e-value threshold
  python3 run_pipeline.py upstream Drosophila_melanogaster -e 1e-10
  
  # Filter existing tblout files
  python3 run_pipeline.py filter Drosophila_melanogaster
  
  # Use custom config file
  python3 run_pipeline.py upstream Drosophila_melanogaster -c custom_config.yaml

The pipeline will:
1. Read configuration from config.yaml (e-value thresholds, excluded genes)
2. Process tblout.gz files in genomes/species_name/
3. Apply filters based on species taxonomy
4. Generate output files in the species directory
        """
    )
    
    parser.add_argument('command', choices=['upstream', 'filter'],
                       help='Pipeline command to execute')
    parser.add_argument('species', 
                       help='Species name (must match directory in genomes/)')
    parser.add_argument('-c', '--config', default='config.yaml',
                       help='Configuration file (default: config.yaml)')
    parser.add_argument('-u', '--upstream', type=int, default=None,
                       help='Upstream distance for extraction (default: from config.yaml)')
    parser.add_argument('-e', '--evalue', type=float, default=1e-7,
                       help='Default e-value cutoff (default: 1e-7)')
    parser.add_argument('-o', '--output',
                       help='Output file (default: auto-generated based on command)')
    parser.add_argument('-i', '--input',
                       help='Input file (for filter command, default: auto-detect)')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Check if species directory exists
    species_dir = Path("genomes") / args.species
    if not species_dir.exists():
        print(f"Error: Species directory '{species_dir}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Execute appropriate command
    if args.command == 'upstream':
        run_upstream_pipeline(
            args.species, 
            config,
            upstream_distance=args.upstream,
            default_evalue=args.evalue,
            output_file=args.output
        )
    elif args.command == 'filter':
        run_filter_pipeline(
            args.species,
            config,
            input_file=args.input,
            output_file=args.output,
            default_evalue=args.evalue
        )

if __name__ == "__main__":
    main()