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
import tempfile
import pandas as pd
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

def get_full_taxonomy(species_name: str) -> Dict[str, str]:
    """Get full taxonomic information for a species from genomes_annotated.tsv"""
    # Clean species name format (convert underscores to spaces for lookup)
    species_lookup = species_name.replace('_', ' ')

    # Try to load the annotated genomes file
    try:
        df = pd.read_csv('data/genomes_annotated.tsv', sep='\t', header=None)
        # Column indices based on the file structure:
        # 0: species, 1: assembly_level, 2: url,
        # 3: phylum, 4: class, 5: order, 6: suborder, 7: infraorder,
        # 8: superfamily, 9: family, 10: genus, 11: subgenus

        # Find the species
        species_row = df[df[0] == species_lookup]

        if not species_row.empty:
            row = species_row.iloc[0]
            taxonomy = {
                'species': species_name,
                'genus': str(row[10]) if pd.notna(row[10]) else None,
                'subgenus': str(row[11]) if pd.notna(row[11]) else None,
                'family': str(row[9]) if pd.notna(row[9]) else None,
                'superfamily': str(row[8]) if pd.notna(row[8]) else None,
                'infraorder': str(row[7]) if pd.notna(row[7]) else None,
                'suborder': str(row[6]) if pd.notna(row[6]) else None,
                'order': str(row[5]) if pd.notna(row[5]) else None,
                'class': str(row[4]) if pd.notna(row[4]) else None,
                'phylum': str(row[3]) if pd.notna(row[3]) else None
            }
            # Remove NA values
            taxonomy = {k: v for k, v in taxonomy.items() if v and v != 'NA'}
            return taxonomy

    except Exception as e:
        print(f"Warning: Could not load taxonomy data: {e}", file=sys.stderr)

    return {'species': species_name}

def get_taxonomy_info(species_name: str) -> Optional[str]:
    """Get taxonomic group for a species (e.g., Insecta for Drosophila)"""
    taxonomy = get_full_taxonomy(species_name)
    # Return class level taxonomy for backward compatibility
    return taxonomy.get('class', 'Insecta')

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
    
    # Check required files - need at least one of tblout.gz or blast.gz files
    required_files = ["genome.fna.gz"]
    if not check_required_files(species_dir, required_files):
        sys.exit(1)
    
    # Check if genome.tsv exists, if not, generate it
    genome_tsv = species_dir / "genome.tsv"
    if not genome_tsv.exists():
        print(f"genome.tsv not found, generating it from genome.fna.gz", file=sys.stderr)
        # Add command to generate genome.tsv at the beginning
        genome_fna = species_dir / "genome.fna.gz"
        generate_tsv_cmd = f"zcat {genome_fna} | ./scripts/get_sequence_lengths.py > {genome_tsv}"
        
        # Execute the command to generate genome.tsv
        result = subprocess.run(generate_tsv_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error generating genome.tsv: {result.stderr}", file=sys.stderr)
            sys.exit(1)
        print(f"Generated {genome_tsv}", file=sys.stderr)
    
    # Check for at least one input file type
    tblout_check = list(species_dir.glob("*tblout.gz"))
    blast_check = list(species_dir.glob("*.blast.gz"))
    if not tblout_check and not blast_check:
        print(f"Error: No tblout.gz or blast.gz files found in {species_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Use upstream distance from config if not specified
    if upstream_distance is None:
        upstream_distance = config.get('upstream_length', 100)
    
    # Get strict mode setting from config
    upstream_strict = config.get('upstream_strict', False)
    
    # Get similarity threshold from config
    similarity_threshold = config.get('similarity_threshold', 99)
    
    # Get low complexity thresholds from config
    low_complexity = config.get('low_complexity', {})
    homopolymer_threshold = low_complexity.get('homopolymer', 0.5)
    dinuc_threshold = low_complexity.get('dinucleotide', 0.75)
    trinuc_threshold = low_complexity.get('trinucleotide', 0.75)
    
    # Build output path
    if output_file is None:
        output_file = species_dir / "upstream.fa"
    
    # Get configuration
    excluded_genes = get_excluded_genes(config, species_name)
    evalue_filters = build_evalue_filters(config)
    
    # Build pipeline command
    pipeline_parts = []
    
    # 1. Extract from tblout and blast files
    tblout_files = list(species_dir.glob("*tblout.gz"))
    blast_files = list(species_dir.glob("*.blast.gz"))
    
    # Create a temporary combined GFF3 file in /tmp
    temp_gff = tempfile.mktemp(suffix='.gff3')
    
    # Build commands to create the temporary GFF3 file
    gff_creation_commands = []
    
    # Add tblout files
    if tblout_files:
        if len(tblout_files) == 1:
            gff_creation_commands.append(f"zcat {tblout_files[0]} | ./scripts/tblout_to_gff3.py >> {temp_gff}")
        else:
            gff_creation_commands.append(f"zcat {' '.join(str(f) for f in tblout_files)} | ./scripts/tblout_to_gff3.py >> {temp_gff}")
    
    # Add blast files with gene name mapping from config
    blast_queries = config.get('blast_queries', {})
    for blast_file in blast_files:
        # Extract query name from filename (e.g., noeCR34335.blast.gz -> noeCR34335)
        query_name = blast_file.stem.replace('.blast', '')
        gene_name = blast_queries.get(query_name, query_name)
        gff_creation_commands.append(f"zcat {blast_file} | ./scripts/blast_to_gff3.py -g {gene_name} >> {temp_gff}")
    
    # Create the combined GFF3 file first, then start the main pipeline
    pipeline_parts.append(f"({'; '.join(gff_creation_commands)}; cat {temp_gff}; rm {temp_gff})")
    
    # 2. No additional GFF3 conversion needed (already done above)
    
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
    extract_cmd = f"./scripts/extract_upstream_regions.py -u {upstream_distance} -l {species_dir}/genome.tsv"
    if upstream_strict:
        extract_cmd += " --strict"
    pipeline_parts.append(extract_cmd)
    
    # 6. Convert to FASTA
    pipeline_parts.append(f"./scripts/gff3_to_fasta.py -g {species_dir}/genome.fna.gz -s {species_name}")
    
    # 7. Filter nearly identical sequences
    pipeline_parts.append(f"./scripts/filter_identical_sequences.py -s {similarity_threshold}")
    
    # 8. Filter low complexity sequences
    filter_low_cmd = f"./scripts/filter_low_complexity.py --homopolymer {homopolymer_threshold} --dinucleotide {dinuc_threshold} --trinucleotide {trinuc_threshold}"
    pipeline_parts.append(filter_low_cmd)

    # 9. Filter non-standard nucleotides (only keep A, C, G, T)
    pipeline_parts.append("./scripts/filter_non_standard_nucleotides.py")

    # 10. Sort sequences
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

def get_promoter_settings(species_name: str, config: Dict) -> Dict:
    """Get promoter annotation settings based on species taxonomy"""
    promoter_config = config.get('promoter_annotation', {})
    default_settings = promoter_config.get('default', {})

    # Get full taxonomy for the species
    taxonomy = get_full_taxonomy(species_name)

    # Check taxonomy levels in order of specificity (species -> genus -> family -> order -> suborder -> infraorder -> class)
    taxonomy_levels = ['species', 'genus', 'family', 'order', 'suborder', 'infraorder', 'class']

    for level in taxonomy_levels:
        taxon = taxonomy.get(level)
        if taxon and taxon in promoter_config:
            settings = promoter_config[taxon].copy()
            # Fill in any missing values with defaults
            for key, value in default_settings.items():
                if key not in settings:
                    settings[key] = value
            return settings

    # Return default settings if no taxonomic match found
    return default_settings

def run_annotate_pipeline(species_name: str, config: Dict, input_file: Optional[str] = None,
                         output_file: Optional[str] = None):
    """Run promoter annotation pipeline using align_promoters_v2.py"""

    species_dir = Path("genomes") / species_name

    # Determine input file
    if input_file is None:
        input_file = species_dir / "upstream.fa"

    # Check if input file exists
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' not found", file=sys.stderr)
        sys.exit(1)

    # Build output path
    if output_file is None:
        output_file = species_dir / "annotated.fa"

    # Check ncRNA diversity before proceeding
    print(f"Checking ncRNA diversity for {species_name}...", file=sys.stderr)
    flag_file = species_dir / "insufficient_diversity.flag"

    diversity_check_cmd = f"python3 ./scripts/check_ncrna_diversity.py {input_file} --write-flag {flag_file}"

    try:
        result = subprocess.run(diversity_check_cmd, shell=True, capture_output=True, text=True)
        # Print the diversity check output
        if result.stderr:
            print(result.stderr, file=sys.stderr)

        if result.returncode != 0:
            print(f"ERROR: Insufficient ncRNA diversity for annotation", file=sys.stderr)
            print(f"  Flag file written to: {flag_file}", file=sys.stderr)
            print(f"  Skipping annotation for {species_name}", file=sys.stderr)
            sys.exit(1)

        # Remove flag file if it exists and we now have sufficient diversity
        if flag_file.exists():
            flag_file.unlink()
            print(f"  Removed old flag file: {flag_file}", file=sys.stderr)

    except Exception as e:
        print(f"Error checking diversity: {e}", file=sys.stderr)
        sys.exit(1)

    # Get promoter settings for this species
    settings = get_promoter_settings(species_name, config)

    if not settings:
        print(f"Error: No promoter annotation settings found for {species_name}", file=sys.stderr)
        sys.exit(1)

    # Build align_promoters_v2.py command - align_promoters_v2.py reads from stdin and writes to stdout
    cmd_parts = ["cat", str(input_file), "|", "./scripts/align_promoters_v2.py"]

    # Add specific kmers if specified (required parameter)
    if 'kmers' in settings:
        kmer_list = ','.join(settings['kmers'])
        cmd_parts.extend(["--kmer", kmer_list])
    else:
        # Default k-mer for insects if not specified
        cmd_parts.extend(["--kmer", "TTCYCAA"])

    # Add promoter length if specified
    if 'promoter_length' in settings:
        cmd_parts.extend(["-l", str(settings['promoter_length'])])

    # Add search range if specified
    if 'search_start' in settings:
        cmd_parts.extend(["-s", str(settings['search_start'])])
    if 'search_end' in settings:
        cmd_parts.extend(["-e", str(settings['search_end'])])

    # Add promoter shift if specified
    if 'promoter_shift' in settings:
        cmd_parts.extend(["-ps", str(settings['promoter_shift'])])

    # Add output redirection
    cmd_parts.extend([">", str(output_file)])

    full_command = " ".join(cmd_parts)

    # Print command for transparency
    print(f"Executing annotate pipeline for {species_name}:", file=sys.stderr)
    print(f"  Command: {full_command}", file=sys.stderr)
    print(f"  Input: {input_file}", file=sys.stderr)
    print(f"  Output: {output_file}", file=sys.stderr)

    # Show applied settings
    taxonomy = get_full_taxonomy(species_name)
    applied_taxon = None
    for level in ['species', 'genus', 'family', 'order', 'suborder', 'infraorder', 'class']:
        taxon = taxonomy.get(level)
        if taxon and taxon in config.get('promoter_annotation', {}):
            applied_taxon = f"{level}: {taxon}"
            break

    if applied_taxon:
        print(f"  Using settings from {applied_taxon}", file=sys.stderr)
    else:
        print(f"  Using default settings", file=sys.stderr)

    # Execute the command
    try:
        result = subprocess.run(full_command, shell=True, capture_output=False, text=True)
        if result.returncode != 0:
            print(f"Error: Annotate pipeline failed with exit code {result.returncode}", file=sys.stderr)
            sys.exit(1)
        print(f"Annotate pipeline completed successfully. Output written to {output_file}", file=sys.stderr)
    except Exception as e:
        print(f"Error executing annotate pipeline: {e}", file=sys.stderr)
        sys.exit(1)

def run_screen_pipeline(species_name: str, config: Dict, cutoff_strategy: str = 'minSum',
                       min_score: float = 0.8, output_file: Optional[str] = None):
    """Run genome screening for Pol II and Pol III promoters using MATCH algorithm"""

    species_dir = Path("genomes") / species_name

    # Check if annotated.fa exists
    annotated_file = species_dir / "annotated.fa"
    if not annotated_file.exists():
        print(f"Error: annotated.fa not found for {species_name}", file=sys.stderr)
        print(f"  Run 'annotate' command first to create aligned sequences", file=sys.stderr)
        sys.exit(1)

    # Check if genome exists
    genome_file = species_dir / "genome.fna.gz"
    if not genome_file.exists():
        print(f"Error: genome.fna.gz not found for {species_name}", file=sys.stderr)
        sys.exit(1)

    # Build output path
    if output_file is None:
        output_file = species_dir / "promoters.gff3"

    print(f"Screening genome for {species_name}...", file=sys.stderr)
    print(f"  Input alignment: {annotated_file}", file=sys.stderr)
    print(f"  Genome: {genome_file}", file=sys.stderr)
    print(f"  Output: {output_file}", file=sys.stderr)

    # Step 1: Build PWMs from annotated sequences
    print("\nStep 1: Building PWMs from aligned sequences...", file=sys.stderr)
    pwm_prefix = species_dir / "pwm"
    build_pwm_cmd = f"./scripts/build_pwm_from_alignment.py {annotated_file} -o {pwm_prefix}"

    try:
        result = subprocess.run(build_pwm_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error building PWMs: {result.stderr}", file=sys.stderr)
            sys.exit(1)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
    except Exception as e:
        print(f"Error running PWM builder: {e}", file=sys.stderr)
        sys.exit(1)

    # Check which PWM files were created
    pol2_pwm_file = f"{pwm_prefix}_pol2.json"
    pol3_pwm_file = f"{pwm_prefix}_pol3.json"

    pol2_exists = Path(pol2_pwm_file).exists()
    pol3_exists = Path(pol3_pwm_file).exists()

    if not pol2_exists and not pol3_exists:
        print("Error: No PWMs could be built from aligned sequences", file=sys.stderr)
        print("  Check that annotated.fa contains properly aligned sequences", file=sys.stderr)
        sys.exit(1)

    # Step 2: Screen genome with PWMs
    print("\nStep 2: Screening genome with fast MATCH algorithm...", file=sys.stderr)
    screen_cmd = f"./scripts/screen_genome_promoters_fast.py {genome_file}"

    if pol2_exists:
        screen_cmd += f" --pol2-pwm {pol2_pwm_file}"
        print(f"  Using Pol II PWM", file=sys.stderr)

    if pol3_exists:
        screen_cmd += f" --pol3-pwm {pol3_pwm_file}"
        print(f"  Using Pol III PWM", file=sys.stderr)

    screen_cmd += f" -o {output_file} -c {cutoff_strategy} -s {min_score}"

    # Add summary output
    summary_file = species_dir / "promoters_summary.txt"
    screen_cmd += f" --summary {summary_file}"

    print(f"  Cutoff strategy: {cutoff_strategy}", file=sys.stderr)
    print(f"  Minimum score: {min_score}", file=sys.stderr)

    try:
        result = subprocess.run(screen_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error screening genome: {result.stderr}", file=sys.stderr)
            sys.exit(1)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
    except Exception as e:
        print(f"Error running genome screening: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"\nScreening completed successfully!", file=sys.stderr)
    print(f"  Predictions: {output_file}", file=sys.stderr)
    print(f"  Summary: {summary_file}", file=sys.stderr)

def run_promoter_pipeline(species_name: str, config: Dict, input_file: Optional[str] = None):
    """Run promoter consensus analysis pipeline"""

    species_dir = Path("genomes") / species_name

    # Determine input file
    if input_file is None:
        input_file = species_dir / "annotated.fa"

    # Check if input file exists
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' not found", file=sys.stderr)
        print(f"  Run 'annotate' command first to create aligned sequences", file=sys.stderr)
        sys.exit(1)

    print(f"Analyzing promoter consensus for {species_name}...", file=sys.stderr)
    print(f"  Input: {input_file}", file=sys.stderr)

    # Run promoter analysis - output goes directly to stdout
    analysis_cmd = f"python3 ./scripts/analyze_promoter_consensus.py {input_file}"

    try:
        result = subprocess.run(analysis_cmd, shell=True, capture_output=False, text=True)
        if result.returncode != 0:
            print(f"Error: Promoter analysis failed with exit code {result.returncode}", file=sys.stderr)
            sys.exit(1)
    except Exception as e:
        print(f"Error running promoter analysis: {e}", file=sys.stderr)
        sys.exit(1)

def run_status_pipeline(config: Dict, show_taxonomy: bool = False, show_rare: bool = False,
                        taxonomy_level: str = 'all', species_filter: Optional[str] = None):
    """Run status analysis showing data collection progress"""

    print("=" * 80, file=sys.stderr)
    print("DATA COLLECTION STATUS", file=sys.stderr)
    print("=" * 80, file=sys.stderr)

    # Run species status collection
    print("\nCollecting species status information...", file=sys.stderr)
    status_cmd = "./scripts/collect_species_status.py -f table"
    if species_filter:
        status_cmd += f" -s {species_filter}"

    try:
        result = subprocess.run(status_cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print("\n=== Species Status Table ===")
            print(result.stdout)
        else:
            print(f"Error collecting species status: {result.stderr}", file=sys.stderr)
    except Exception as e:
        print(f"Error running species status: {e}", file=sys.stderr)

    # Run summary statistics
    print("\nGenerating summary statistics...", file=sys.stderr)
    summary_cmd = "./scripts/collect_species_status.py -f summary"

    try:
        result = subprocess.run(summary_cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print("\n=== Summary Statistics ===")
            print(result.stdout)
    except Exception:
        pass

    # Run taxonomic analysis if requested
    if show_taxonomy:
        print("\nAnalyzing ncRNA distribution by taxonomy...", file=sys.stderr)
        tax_cmd = f"./scripts/analyze_ncrna_taxonomy.py -l {taxonomy_level}"
        if show_rare:
            tax_cmd += " --rare"

        try:
            result = subprocess.run(tax_cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                print("\n=== Taxonomic Analysis ===")
                print(result.stdout)
            else:
                print(f"Error in taxonomic analysis: {result.stderr}", file=sys.stderr)
        except Exception as e:
            print(f"Error running taxonomic analysis: {e}", file=sys.stderr)

    print("\n" + "=" * 80, file=sys.stderr)
    print("Use 'python3 run_pipeline.py status --taxonomy' for taxonomic analysis", file=sys.stderr)
    print("Use 'python3 run_pipeline.py status --rare' to show rare ncRNAs", file=sys.stderr)
    print("=" * 80, file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Run RNA analysis pipelines using configuration from config.yaml',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Commands:
  upstream    Extract upstream regions and generate FASTA file
  filter      Apply filtering to tblout or GFF3 files
  annotate    Annotate promoter regions using taxonomic-specific settings
  promoter    Analyze promoter consensus and Pol II/III differences
  screen      Screen genome for promoters using MATCH algorithm
  status      Show data collection status and ncRNA distribution

Combined Commands (use + to chain):
  upstream+annotate    Extract upstream regions then annotate promoters
  upstream+annotate+promoter   Extract, annotate, and analyze promoter consensus
  upstream+annotate+screen    Complete pipeline from extraction to screening

Examples:
  # Extract upstream regions for Drosophila melanogaster
  python3 run_pipeline.py upstream Drosophila_melanogaster

  # Extract and annotate in one command
  python3 run_pipeline.py upstream+annotate Drosophila_melanogaster

  # Complete pipeline from extraction to screening
  python3 run_pipeline.py upstream+annotate+screen Drosophila_melanogaster

  # Extract with custom upstream distance
  python3 run_pipeline.py upstream Drosophila_melanogaster -u 200

  # Extract with custom e-value threshold
  python3 run_pipeline.py upstream Drosophila_melanogaster -e 1e-10

  # Filter existing tblout files
  python3 run_pipeline.py filter Drosophila_melanogaster

  # Annotate promoter regions (uses upstream.fa as input)
  python3 run_pipeline.py annotate Drosophila_melanogaster

  # Analyze promoter consensus (uses annotated.fa as input)
  python3 run_pipeline.py promoter Drosophila_melanogaster

  # Screen genome for promoters (requires annotated.fa)
  python3 run_pipeline.py screen Drosophila_melanogaster

  # Show status for all species
  python3 run_pipeline.py status

  # Show status with taxonomic analysis
  python3 run_pipeline.py status --taxonomy

  # Show status for specific taxonomic level
  python3 run_pipeline.py status --taxonomy --level family

  # Use custom config file
  python3 run_pipeline.py upstream Drosophila_melanogaster -c custom_config.yaml

The pipeline will:
1. Read configuration from config.yaml (e-value thresholds, excluded genes, promoter settings)
2. Process tblout.gz files in genomes/species_name/
3. Apply filters based on species taxonomy
4. Generate output files in the species directory
        """
    )

    parser.add_argument('command',
                       help='Pipeline command to execute (single: upstream, filter, annotate, screen, status; or combined: upstream+annotate, upstream+annotate+screen)')
    parser.add_argument('species', nargs='?',
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
    parser.add_argument('--taxonomy', action='store_true',
                       help='Show taxonomic analysis (for status command)')
    parser.add_argument('--rare', action='store_true',
                       help='Show rare ncRNAs (for status command)')
    parser.add_argument('--level', choices=['genus', 'family', 'order', 'class', 'all'],
                       default='all',
                       help='Taxonomic level for analysis (for status command)')
    parser.add_argument('--cutoff-strategy', choices=['minFN', 'minFP', 'minSum'],
                       default='minSum',
                       help='MATCH cutoff strategy (for screen command, default: minSum)')
    parser.add_argument('--min-score', type=float, default=0.8,
                       help='Minimum score for promoter predictions (for screen command, default: 0.8)')

    args = parser.parse_args()

    # Parse command - handle combined commands with '+'
    commands = args.command.split('+')
    valid_commands = ['upstream', 'filter', 'annotate', 'promoter', 'screen', 'status']

    # Validate all commands
    for cmd in commands:
        if cmd not in valid_commands:
            print(f"Error: Invalid command '{cmd}'. Valid commands: {', '.join(valid_commands)}", file=sys.stderr)
            sys.exit(1)

    # Load configuration
    config = load_config(args.config)

    # For status command, species is optional
    if 'status' not in commands and not args.species:
        print("Error: Species name is required for this command", file=sys.stderr)
        sys.exit(1)

    # Check if species directory exists (except for status command)
    if 'status' not in commands:
        species_dir = Path("genomes") / args.species
        if not species_dir.exists():
            print(f"Error: Species directory '{species_dir}' not found", file=sys.stderr)
            sys.exit(1)

    # Execute commands in sequence
    print(f"Executing pipeline: {' -> '.join(commands)}", file=sys.stderr)

    for i, command in enumerate(commands):
        if len(commands) > 1:
            print(f"\n{'='*60}", file=sys.stderr)
            print(f"STEP {i+1}/{len(commands)}: {command.upper()}", file=sys.stderr)
            print(f"{'='*60}", file=sys.stderr)

        if command == 'upstream':
            run_upstream_pipeline(
                args.species,
                config,
                upstream_distance=args.upstream,
                default_evalue=args.evalue,
                output_file=args.output
            )
        elif command == 'filter':
            run_filter_pipeline(
                args.species,
                config,
                input_file=args.input,
                output_file=args.output,
                default_evalue=args.evalue
            )
        elif command == 'annotate':
            run_annotate_pipeline(
                args.species,
                config,
                input_file=args.input,
                output_file=args.output
            )
        elif command == 'promoter':
            run_promoter_pipeline(
                args.species,
                config,
                input_file=args.input
            )
        elif command == 'screen':
            run_screen_pipeline(
                args.species,
                config,
                cutoff_strategy=args.cutoff_strategy,
                min_score=args.min_score,
                output_file=args.output
            )
        elif command == 'status':
            run_status_pipeline(
                config,
                show_taxonomy=args.taxonomy,
                show_rare=args.rare,
                taxonomy_level=args.level,
                species_filter=args.species
            )

    if len(commands) > 1:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"PIPELINE COMPLETE: All {len(commands)} steps finished successfully!", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)

if __name__ == "__main__":
    main()