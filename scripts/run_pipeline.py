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

    # Add promoter start if specified
    if 'promoter_start' in settings:
        cmd_parts.extend(["-ps", str(settings['promoter_start'])])

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

def run_promoterwindow_pipeline(species_name: str, config: Dict, input_file: Optional[str] = None, eval_mode: bool = False):
    """Run optimal promoter window analysis pipeline"""
    species_dir = Path("genomes") / species_name

    # Step 1: Generate core motif alignment with fixed window
    print("Step 1: Generating core motif alignment...", file=sys.stderr)

    # Get promoter settings for this species
    settings = get_promoter_settings(species_name, config)

    if not settings or 'kmers' not in settings or not settings['kmers']:
        print(f"Error: No k-mers found for {species_name}", file=sys.stderr)
        sys.exit(1)

    # Use the first k-mer from the settings
    kmer = settings['kmers'][0]
    kmer_length = len(kmer)
    print(f"Using k-mer: {kmer} (length: {kmer_length})", file=sys.stderr)

    # Run annotate pipeline with regular parameters to get proper promoter windows
    temp_annotated = species_dir / "annotated_for_window.fa"

    # Use all settings from the promoter configuration
    annotate_cmd = [
        "python3", "scripts/align_promoters_v2.py",
        "--kmer", kmer,
        "-l", str(settings.get('promoter_length', 20)),  # Use configured promoter length
        "-ps", str(settings.get('promoter_start', 3)),   # Use configured promoter start
        "-s", "-40", "-e", "-60"
        # Removed -v for verbose output
    ]

    print(f"Running: {' '.join(annotate_cmd)} < {species_dir / 'upstream.fa'}", file=sys.stderr)
    print(f"  Using promoter length: {settings.get('promoter_length', 20)}, start: {settings.get('promoter_start', 3)}", file=sys.stderr)

    with open(species_dir / "upstream.fa", 'r') as input_f:
        with open(temp_annotated, 'w') as output_f:
            with open(os.devnull, 'w') as devnull:
                result = subprocess.run(annotate_cmd, stdin=input_f, stdout=output_f, stderr=devnull, text=True)

    if result.returncode != 0:
        print(f"Error in alignment step", file=sys.stderr)
        sys.exit(1)

    # Step 2: Analyze optimal window
    if eval_mode:
        print("Step 2: Evaluating multiple thresholds...", file=sys.stderr)
    else:
        print("Step 2: Analyzing optimal promoter window...", file=sys.stderr)

    output_file = species_dir / "promoter_window_analysis.yaml"

    window_cmd = [
        "python3", "scripts/find_optimal_promoter_window.py",
        str(temp_annotated),
        "--kmer", kmer,
        "-c", "config.yaml"
    ]

    if eval_mode:
        window_cmd.append("--eval")
    else:
        window_cmd.append("--verbose")

    print(f"Running: {' '.join(window_cmd)}", file=sys.stderr)
    result = subprocess.run(window_cmd)

    if result.returncode != 0:
        print("Error in window analysis", file=sys.stderr)
        sys.exit(1)

    # Clean up temporary file
    temp_annotated.unlink(missing_ok=True)
    # print(f"Temporary alignment file kept at: {temp_annotated}", file=sys.stderr)  # For debugging

    print(f"Promoter window analysis completed.", file=sys.stderr)

def run_noeCR34335_pipeline(species_name: str, config: Dict, output_file: Optional[str] = None):
    """Extract noeCR34335 lncRNA hits with promoter annotations for publication"""

    species_dir = Path("genomes") / species_name

    # Check if annotated.fa exists
    annotated_file = species_dir / "annotated.fa"
    if not annotated_file.exists():
        print(f"Error: annotated.fa not found for {species_name}", file=sys.stderr)
        print(f"  Run 'annotate' command first to create aligned sequences", file=sys.stderr)
        sys.exit(1)

    # Check if noeCR34335.blast.gz exists
    blast_file = species_dir / "noeCR34335.blast.gz"
    if not blast_file.exists():
        print(f"Error: noeCR34335.blast.gz not found for {species_name}", file=sys.stderr)
        sys.exit(1)

    # Check if genome exists
    genome_file = species_dir / "genome.fna.gz"
    if not genome_file.exists():
        print(f"Error: genome.fna.gz not found for {species_name}", file=sys.stderr)
        sys.exit(1)

    # Build output path
    results_dir = Path("results") / species_name
    results_dir.mkdir(parents=True, exist_ok=True)

    if output_file is None:
        output_file = results_dir / "lncrna.fa"

    print(f"Extracting noeCR34335 lncRNA data for {species_name}...", file=sys.stderr)
    print(f"  Promoter info: {annotated_file}", file=sys.stderr)
    print(f"  BLAST hits: {blast_file}", file=sys.stderr)
    print(f"  Genome: {genome_file}", file=sys.stderr)
    print(f"  Output FASTA: {output_file}", file=sys.stderr)

    # Step 1: Load species abbreviations
    import re
    from Bio import SeqIO
    import pandas as pd

    # Load species abbreviations
    abbreviations = {}
    try:
        df = pd.read_csv('data/species_abbreviations.tsv', sep='\t', header=None)
        for _, row in df.iterrows():
            species_underscore = row[1]
            abbrev = row[2]
            abbreviations[species_underscore] = abbrev
        print(f"  Loaded {len(abbreviations)} species abbreviations", file=sys.stderr)
    except Exception as e:
        print(f"  Warning: Could not load species abbreviations: {e}", file=sys.stderr)

    # Get abbreviation for current species
    species_abbrev = abbreviations.get(species_name, species_name)
    print(f"  Using abbreviation: {species_abbrev}", file=sys.stderr)

    hits_with_promoters = set()
    upstream_with_promoters = []
    noecr_counter = 1

    # Parse the annotated.fa file to find all sequences with promoters
    for record in SeqIO.parse(annotated_file, "fasta"):
        # Check if ANY sequence has a promoter (capital letters in the sequence)
        seq_str = str(record.seq)
        promoter_match = re.search(r'[ACGT]{10,}', seq_str)
        if promoter_match:
            # This sequence has a promoter, add it to upstream collection
            upstream_with_promoters.append(record)

            # If this is also a noeCR34335 hit, track it separately for lncRNA extraction
            if "noeCR34335" in record.id:
                parts = record.id.split('|')
                if len(parts) >= 6:
                    chrom = parts[2]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[5]
                    hits_with_promoters.add((chrom, start, end, strand, noecr_counter))
                    print(f"    Found noeCR34335 promoter for upstream region {chrom}:{start}-{end}({strand}) -> {species_abbrev}{noecr_counter}", file=sys.stderr)
                    noecr_counter += 1

    if not hits_with_promoters:
        print(f"Warning: No noeCR34335 hits with promoters found in {annotated_file}", file=sys.stderr)
        return

    print(f"  Found {len(upstream_with_promoters)} total sequences with promoters", file=sys.stderr)
    print(f"  Found {len(hits_with_promoters)} noeCR34335 hits with promoters", file=sys.stderr)

    # Step 2: Convert BLAST file to GFF3 and filter for hits with promoters
    temp_gff = tempfile.mktemp(suffix='.gff3')
    filtered_gff = tempfile.mktemp(suffix='.gff3')

    # Convert BLAST to GFF3
    print(f"  Converting BLAST hits to GFF3...", file=sys.stderr)
    blast_to_gff_cmd = f"zcat {blast_file} | ./scripts/blast_to_gff3.py -g noeCR34335 > {temp_gff}"

    result = subprocess.run(blast_to_gff_cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error converting BLAST to GFF3: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    # Step 3: Filter GFF3 for only those hits that have promoters
    print(f"  Filtering for ncRNAs with promoters...", file=sys.stderr)
    matched_count = 0

    with open(temp_gff, 'r') as infile, open(filtered_gff, 'w') as outfile:
        outfile.write("##gff-version 3\n")
        outfile.write(f"##source-version noeCR34335 pipeline 1.0\n")
        outfile.write(f"##species {species_name}\n")

        for line in infile:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) >= 9:
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]

                # Check if this ncRNA has a corresponding upstream promoter
                # The upstream region should be adjacent to this ncRNA
                # For - strand: upstream starts at ncRNA_end+1 and ends at ncRNA_end+100
                # For + strand: upstream ends at ncRNA_start-1 and starts at ncRNA_start-100
                has_promoter = False
                noecr_id = ""
                for (up_chrom, up_start, up_end, up_strand, counter) in hits_with_promoters:
                    if chrom == up_chrom and strand == up_strand:
                        # Check if this upstream region is adjacent to this ncRNA
                        # Note: BLAST coordinates might be swapped for minus strand
                        actual_start = min(start, end)
                        actual_end = max(start, end)

                        if strand == '-' and up_start == actual_end + 1 and up_end == actual_end + 100:
                            has_promoter = True
                            noecr_id = f"{species_abbrev}{counter}"
                            print(f"      Matched - strand: ncRNA {actual_start}-{actual_end}, upstream {up_start}-{up_end} -> {noecr_id}", file=sys.stderr)
                            break
                        elif strand == '+' and up_end == actual_start - 1 and up_start == actual_start - 100:
                            has_promoter = True
                            noecr_id = f"{species_abbrev}{counter}"
                            print(f"      Matched + strand: ncRNA {actual_start}-{actual_end}, upstream {up_start}-{up_end} -> {noecr_id}", file=sys.stderr)
                            break

                if has_promoter:
                    # Replace the line with new ID and store mapping
                    parts = line.strip().split('\t')
                    parts[8] = f"ID={noecr_id};Name=noeCR34335"
                    outfile.write('\t'.join(parts) + '\n')

                    # Store coordinate mapping for ID replacement
                    coord_key = f"{chrom}:{start}-{end}({strand})"
                    if 'blast_to_id_mapping' not in locals():
                        blast_to_id_mapping = {}
                    blast_to_id_mapping[coord_key] = noecr_id

                    matched_count += 1

    print(f"  Matched {matched_count} ncRNAs with promoters", file=sys.stderr)

    # Step 4: Extract the ncRNA sequences from the genome with poly-T extension
    print(f"  Extracting ncRNA sequences from genome with poly-T extension...", file=sys.stderr)

    temp_output = tempfile.mktemp(suffix='.fa')
    incomplete_file = results_dir / 'lncrna.incomplete'
    extract_cmd = f"./scripts/extract_with_polyt_extension.py -g {genome_file} -s {species_name} --incomplete {incomplete_file} < {filtered_gff} > {temp_output}"

    result = subprocess.run(extract_cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error extracting sequences: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    # Use the blast_to_id_mapping created during filtering
    if 'blast_to_id_mapping' not in locals():
        blast_to_id_mapping = {}

    print(f"  Using {len(blast_to_id_mapping)} coordinate mappings for ID replacement", file=sys.stderr)

    with open(temp_output, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Parse header to extract coordinates
                # Format: >species-number|gene|chrom|start|end|strand
                header = line.strip()
                if '|' in header:
                    parts = header[1:].split('|')
                    if len(parts) >= 6:
                        chrom = parts[2]
                        start = parts[3]
                        end = parts[4]
                        strand = parts[5]
                        coord_key = f"{chrom}:{start}-{end}({strand})"
                        if coord_key in blast_to_id_mapping:
                            new_id = blast_to_id_mapping[coord_key]
                            outfile.write(f">{new_id}|noeCR34335|{chrom}|{start}|{end}|{strand}\n")
                        else:
                            outfile.write(line)
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)
            else:
                outfile.write(line)

    # Clean up temp file
    try:
        os.unlink(temp_output)
    except:
        pass

    print(f"  Wrote {matched_count} ncRNA sequences to {output_file}", file=sys.stderr)

    # Update incomplete file with mapped IDs if it exists
    if incomplete_file.exists():
        print(f"  Updating incomplete sequence IDs...", file=sys.stderr)
        updated_incomplete = []
        with open(incomplete_file, 'r') as f:
            for line in f:
                original_id = line.strip()
                if original_id:
                    # Extract the sequence number from the original ID (e.g., "Drosophila_burlai-1" -> "1")
                    if '-' in original_id:
                        seq_num = original_id.split('-')[-1]
                        # Find the corresponding new ID
                        new_id = f"{species_abbrev}{seq_num}"
                        updated_incomplete.append(new_id)
                    else:
                        updated_incomplete.append(original_id)

        # Write updated incomplete file
        with open(incomplete_file, 'w') as f:
            for seq_id in updated_incomplete:
                f.write(f"{seq_id}\n")
        print(f"  Updated {len(updated_incomplete)} incomplete sequence IDs", file=sys.stderr)

    # Step 5: Copy filtered upstream sequences to results directory with updated IDs
    upstream_output = results_dir / "upstream.fa"
    print(f"  Copying filtered upstream sequences to {upstream_output}...", file=sys.stderr)

    # Create mapping from coordinates to new IDs
    coord_to_id = {}
    for (chrom, start, end, strand, counter) in hits_with_promoters:
        coord_to_id[(chrom, start, end, strand)] = f"{species_abbrev}{counter}"

    # Write sequences in 2-line format to match original
    with open(upstream_output, 'w') as f:
        for record in upstream_with_promoters:
            # Check if this is a noeCR34335 hit and update ID
            if "noeCR34335" in record.id:
                parts = record.id.split('|')
                if len(parts) >= 6:
                    chrom = parts[2]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[5]
                    new_id = coord_to_id.get((chrom, start, end, strand), None)
                    if new_id:
                        # Create new header with updated ID but keeping noeCR34335 as gene type
                        f.write(f">{new_id}|noeCR34335|{chrom}|{start}|{end}|{strand}\n")
                    else:
                        f.write(f">{record.id}\n")
                else:
                    f.write(f">{record.id}\n")
            else:
                f.write(f">{record.id}\n")
            f.write(f"{str(record.seq)}\n")

    print(f"  Wrote {len(upstream_with_promoters)} upstream sequences with promoters to {upstream_output}", file=sys.stderr)

    # Step 6: Build Pol II and Pol III PWMs and score noeCR34335 promoters
    print(f"\n  Building PWMs for scoring...", file=sys.stderr)

    # Load polymerase classifications from config
    import yaml
    try:
        with open('config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        pol2_types = config['rna_polymerases']['pol2']['types']
        pol3_types = config['rna_polymerases']['pol3']['types']
        print(f"  Pol II types: {pol2_types}", file=sys.stderr)
        print(f"  Pol III types: {pol3_types}", file=sys.stderr)
    except:
        print(f"  Warning: Could not load config.yaml, using default classifications", file=sys.stderr)
        pol2_types = ['U1', 'U2', 'U3', 'U4', 'U4atac', 'U5', 'U7', 'U11', 'U12', 'OrCD1']
        pol3_types = ['Arthropod_7SK', 'RNase_MRP', 'RNaseP_nuc', 'U6', 'U6atac', 'tRNAsec', 'snoRNAMe18SA1806']

    # Create separate files for Pol II and Pol III sequences
    pol2_file = tempfile.mktemp(suffix='_pol2.fa')
    pol3_file = tempfile.mktemp(suffix='_pol3.fa')
    noecr_promoters = []

    pol2_count = 0
    pol3_count = 0
    pol2_sequences = []
    pol3_sequences = []

    with open(pol2_file, 'w') as pol2_f, open(pol3_file, 'w') as pol3_f:
        for record in upstream_with_promoters:
            if "noeCR34335" in record.id:
                # Extract noeCR34335 promoter sequence for scoring
                seq_str = str(record.seq)
                promoter_match = re.search(r'[ACGT]{15,25}', seq_str)
                if promoter_match:
                    promoter_seq = promoter_match.group(0)
                    # Get the new ID for this sequence
                    parts = record.id.split('|')
                    if len(parts) >= 6:
                        chrom = parts[2]
                        start = int(parts[3])
                        end = int(parts[4])
                        strand = parts[5]
                        new_id = coord_to_id.get((chrom, start, end, strand), parts[0])
                        noecr_promoters.append((new_id, promoter_seq))
            else:
                # Classify into Pol II or Pol III based on gene type (excluding noeCR34335)
                gene_type = record.id.split('|')[1] if '|' in record.id else record.id.split('-')[1] if '-' in record.id else "unknown"

                # Extract uppercase nucleotides (the actual promoter sequence)
                promoter_seq = ''.join([c for c in str(record.seq) if c.isupper()])

                if gene_type in pol2_types:
                    pol2_f.write(f">{record.id}\n")
                    pol2_f.write(f"{str(record.seq)}\n")
                    pol2_sequences.append((gene_type, promoter_seq))
                    pol2_count += 1
                elif gene_type in pol3_types and gene_type != 'noeCR34335':
                    pol3_f.write(f">{record.id}\n")
                    pol3_f.write(f"{str(record.seq)}\n")
                    pol3_sequences.append((gene_type, promoter_seq))
                    pol3_count += 1

    print(f"  Found {len(noecr_promoters)} noeCR34335 sequences to score:", file=sys.stderr)
    for seq_id, seq in noecr_promoters:
        print(f"    {seq_id}: {seq}", file=sys.stderr)

    print(f"  Found {pol2_count} Pol II sequences and {pol3_count} Pol III sequences for PWM building", file=sys.stderr)

    if pol2_sequences:
        print(f"  Pol II PWM sequences ({len(pol2_sequences)} sequences):", file=sys.stderr)
        for gene_type, seq in pol2_sequences:
            print(f"    {gene_type}: {seq}", file=sys.stderr)

    if pol3_sequences:
        print(f"  Pol III PWM sequences ({len(pol3_sequences)} sequences, excluding noeCR34335):", file=sys.stderr)
        for gene_type, seq in pol3_sequences:
            print(f"    {gene_type}: {seq}", file=sys.stderr)

    # Build PWMs if we have enough sequences
    pwm_prefix = results_dir / "background_pwm"
    pwms_built = True

    if pol2_count >= 2:
        print(f"  Building Pol II PWM from {pol2_count} sequences...", file=sys.stderr)
        build_pol2_cmd = f"./scripts/build_pwm_from_alignment.py {pol2_file} -o {pwm_prefix}_pol2"
        result = subprocess.run(build_pol2_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  Warning: Pol II PWM building failed: {result.stderr}", file=sys.stderr)
            pwms_built = False
    else:
        print(f"  Warning: Insufficient Pol II sequences ({pol2_count}) for PWM building", file=sys.stderr)
        pwms_built = False

    if pol3_count >= 2:
        print(f"  Building Pol III PWM from {pol3_count} sequences (excluding noeCR34335)...", file=sys.stderr)
        build_pol3_cmd = f"./scripts/build_pwm_from_alignment.py {pol3_file} -o {pwm_prefix}_pol3"
        result = subprocess.run(build_pol3_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  Warning: Pol III PWM building failed: {result.stderr}", file=sys.stderr)
            pwms_built = False
    else:
        print(f"  Warning: Insufficient Pol III sequences ({pol3_count}) for PWM building", file=sys.stderr)
        pwms_built = False

    if pwms_built and noecr_promoters:
        print(f"  Built background PWMs", file=sys.stderr)

        # Score noeCR34335 promoters
        score_output = results_dir / "upstream.score"
        print(f"  Scoring noeCR34335 promoters...", file=sys.stderr)

        with open(score_output, 'w') as score_f:
            # Write header
            score_f.write("ID\tSequence\tPol2_Score\tPol3_Score\n")

            pol2_pwm = f"{pwm_prefix}_pol2_pol2.json"
            pol3_pwm = f"{pwm_prefix}_pol3_pol3.json"

            for seq_id, promoter_seq in noecr_promoters:
                pol2_score = pol3_score = "NA"

                # Score against Pol II PWM
                if os.path.exists(pol2_pwm):
                    pol2_cmd = f"echo '>{seq_id}\\n{promoter_seq}' | python3 ./scripts/score_sequences_with_pwm.py {pol2_pwm}"
                    pol2_result = subprocess.run(pol2_cmd, shell=True, capture_output=True, text=True)
                    if pol2_result.returncode == 0 and pol2_result.stdout.strip():
                        lines = pol2_result.stdout.strip().split('\n')
                        for line in lines:
                            if '\t' in line:
                                parts = line.split('\t')
                                if len(parts) >= 2:
                                    try:
                                        pol2_score = f"{float(parts[1]):.2f}"
                                        break
                                    except:
                                        pass

                # Score against Pol III PWM
                if os.path.exists(pol3_pwm):
                    pol3_cmd = f"echo '>{seq_id}\\n{promoter_seq}' | python3 ./scripts/score_sequences_with_pwm.py {pol3_pwm}"
                    pol3_result = subprocess.run(pol3_cmd, shell=True, capture_output=True, text=True)
                    if pol3_result.returncode == 0 and pol3_result.stdout.strip():
                        lines = pol3_result.stdout.strip().split('\n')
                        for line in lines:
                            if '\t' in line:
                                parts = line.split('\t')
                                if len(parts) >= 2:
                                    try:
                                        pol3_score = f"{float(parts[1]):.2f}"
                                        break
                                    except:
                                        pass

                score_f.write(f"{seq_id}\t{promoter_seq}\t{pol2_score}\t{pol3_score}\n")

            print(f"  Wrote promoter scores to {score_output}", file=sys.stderr)
    else:
        print(f"  Skipping scoring due to PWM building issues", file=sys.stderr)

    # Clean up temporary files
    try:
        os.unlink(temp_gff)
        os.unlink(filtered_gff)
        os.unlink(pol2_file)
        os.unlink(pol3_file)
    except:
        pass

    # Step 7: Score sequences with cmsearch
    print(f"\n  Scoring sequences with cmsearch...", file=sys.stderr)
    noe_scores_file = results_dir / "noe.scores"

    cmsearch_cmd = f"./scripts/score_sequences_with_cmsearch.py {output_file} -o {noe_scores_file} -v"
    result = subprocess.run(cmsearch_cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  Warning: cmsearch scoring failed: {result.stderr}", file=sys.stderr)
    else:
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        print(f"  cmsearch scores written to {noe_scores_file}", file=sys.stderr)

    # Step 8: Generate HTML and PDF reports
    print(f"\n  Generating HTML and PDF reports...", file=sys.stderr)

    html_file = results_dir / "lncrna.html"
    pdf_file = results_dir / "lncrna.pdf"

    # Remove existing files
    for file_path in [html_file, pdf_file]:
        if file_path.exists():
            file_path.unlink()

    # Generate HTML report
    with open(html_file, 'w') as f:
        # HTML header with styles
        f.write("""<!DOCTYPE html>
<html>
<head>
<style>
@page { size: A4 portrait; margin: 1em; }
body { font-family: Arial, Helvetica, sans-serif; }
.alignment { line-height: 1.0; margin: 0; padding: 0; }
.label-green {
    display: inline-block; padding: 1px 1px; background-color: #4CAF50;
    color: white; font-size: 8px; border-radius: 2px;
    font-family: Arial, sans-serif; box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    font-weight: bold;
}
.label-yellow {
    display: inline-block; padding: 1px 1px; background-color: #FFA500;
    color: white; font-size: 8px; border-radius: 2px;
    font-family: Arial, sans-serif; box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    font-weight: bold;
}
.label-red {
    display: inline-block; padding: 1px 1px; background-color: red;
    color: white; font-size: 8px; border-radius: 2px;
    font-family: Arial, sans-serif; box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    font-weight: bold;
}
</style>
</head>
<body>
""")

        # Species information
        species_display = species_name.replace('_', ' ')
        f.write(f"<h2>{species_display}</h2>\n")

        # lncRNA sequences section
        f.write("<h4>Identified lncRNAs</h4>\n")
        f.write('<pre class="alignment" style="font-size:8px;font-family:monospace;line-height:1.0;margin:0;padding:0;">\n')

        # Format lncRNA sequences with scores
        if output_file.exists():
            format_cmd = f"cat {output_file} | python3 ./scripts/html_format_seq.py {score_output} {incomplete_file} {noe_scores_file}"
            result = subprocess.run(format_cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                f.write(result.stdout)
            else:
                f.write("Error formatting lncRNA sequences\n")

        f.write("</pre>\n")

        # Upstream promoter alignment section
        f.write("<h4>Upstream promoter alignment</h4>\n")
        f.write('<pre class="alignment" style="font-size:8px;font-family:monospace;line-height:1.0;margin:0;padding:0;">\n')

        # Format ALL upstream sequences in ClustalW-like alignment format
        if upstream_output.exists():
            # Use the new colorize_alignment script for ALL sequences
            alignment_cmd = f"cat {upstream_output} | python3 ./scripts/colorize_alignment.py"
            result = subprocess.run(alignment_cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                f.write(result.stdout)
            else:
                f.write("Error formatting upstream alignment\n")

        f.write("</pre>\n")

        f.write("</body>\n</html>\n")

    print(f"  Generated HTML report: {html_file}", file=sys.stderr)

    # Generate PDF using headless Chrome
    pdf_cmd = f"google-chrome --headless --disable-gpu --no-pdf-header-footer --print-to-pdf={pdf_file} {html_file}"
    result = subprocess.run(pdf_cmd, shell=True, capture_output=True, text=True)

    if result.returncode == 0 and pdf_file.exists():
        print(f"  Generated PDF report: {pdf_file}", file=sys.stderr)
    else:
        print(f"  Warning: PDF generation failed. HTML report available at {html_file}", file=sys.stderr)

    print(f"\nnoeCR34335 lncRNA extraction completed for {species_name}", file=sys.stderr)

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
  promoterwindow  Find optimal promoter window length based on conservation
  screen      Screen genome for promoters using MATCH algorithm
  status      Show data collection status and ncRNA distribution
  noeCR34335  Extract noeCR34335 lncRNA hits with promoter annotations

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

  # Extract noeCR34335 lncRNA hits with promoters
  python3 run_pipeline.py noeCR34335 Drosophila_melanogaster

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
    parser.add_argument('--eval', action='store_true',
                       help='Evaluate multiple thresholds for promoterwindow command')

    args = parser.parse_args()

    # Parse command - handle combined commands with '+'
    commands = args.command.split('+')
    valid_commands = ['upstream', 'filter', 'annotate', 'promoter', 'promoterwindow', 'screen', 'status', 'noeCR34335']

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
        elif command == 'promoterwindow':
            run_promoterwindow_pipeline(
                args.species,
                config,
                input_file=args.input,
                eval_mode=args.eval
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
        elif command == 'noeCR34335':
            run_noeCR34335_pipeline(
                args.species,
                config,
                output_file=args.output
            )

    if len(commands) > 1:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"PIPELINE COMPLETE: All {len(commands)} steps finished successfully!", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)

if __name__ == "__main__":
    main()