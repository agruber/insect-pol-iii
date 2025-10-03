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
    
    # Check if genome.tsv exists in cache, if not, generate it
    cache_dir = Path("cache")
    cache_dir.mkdir(exist_ok=True)
    genome_tsv = cache_dir / f"{species_name}.tsv"

    if not genome_tsv.exists():
        print(f"Genome lengths not cached, generating {genome_tsv} from genome.fna.gz", file=sys.stderr)
        # Add command to generate genome.tsv in cache
        genome_fna = species_dir / "genome.fna.gz"
        generate_tsv_cmd = f"zcat {genome_fna} | ./scripts/get_sequence_lengths.py > {genome_tsv}"

        # Execute the command to generate genome.tsv
        result = subprocess.run(generate_tsv_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error generating genome lengths: {result.stderr}", file=sys.stderr)
            sys.exit(1)
        print(f"Generated genome lengths cache: {genome_tsv}", file=sys.stderr)
    else:
        print(f"Using cached genome lengths: {genome_tsv}", file=sys.stderr)
    
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
    extract_cmd = f"./scripts/extract_upstream_regions.py -u {upstream_distance} -l {genome_tsv}"
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

    # Check taxonomy levels in order of specificity (species -> genus -> family -> superfamily -> order -> suborder -> infraorder -> class)
    taxonomy_levels = ['species', 'genus', 'family', 'superfamily', 'order', 'suborder', 'infraorder', 'class']

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
    # Handle both old format (direct values) and new format (nested with defaults and RNA-specific)
    if 'search_start' in settings:
        if isinstance(settings['search_start'], dict):
            # New nested format
            if 'default' in settings['search_start']:
                cmd_parts.extend(["-s", str(settings['search_start']['default'])])
            # Add RNA-specific parameters
            for rna_type, value in settings['search_start'].items():
                if rna_type != 'default':
                    cmd_parts.extend(["--rna-params", f"{rna_type}:search_start:{value}"])
        else:
            # Old format - single value
            cmd_parts.extend(["-s", str(settings['search_start'])])

    if 'search_end' in settings:
        if isinstance(settings['search_end'], dict):
            # New nested format
            if 'default' in settings['search_end']:
                cmd_parts.extend(["-e", str(settings['search_end']['default'])])
            # Add RNA-specific parameters
            for rna_type, value in settings['search_end'].items():
                if rna_type != 'default':
                    cmd_parts.extend(["--rna-params", f"{rna_type}:search_end:{value}"])
        else:
            # Old format - single value
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
    for level in ['species', 'genus', 'family', 'superfamily', 'order', 'suborder', 'infraorder', 'class']:
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

def get_assembly_name(species_name: str) -> str:
    """Get assembly name from genome.fna.gz symlink"""
    species_dir = Path("genomes") / species_name
    genome_file = species_dir / "genome.fna.gz"

    try:
        if genome_file.is_symlink():
            # Get the target of the symlink
            target = genome_file.resolve().name
            # Remove .fna.gz extension to get assembly name
            assembly_name = target.replace('.fna.gz', '')
            return assembly_name
        else:
            # If not a symlink, try to extract from filename
            return genome_file.name.replace('.fna.gz', '')
    except:
        return "Unknown assembly"

def trim_polyt_proper(sequence):
    """
    Properly trim polyT from 3' end, including cases like ...ATCGCATTTTTT.

    Steps:
    1. Trim trailing Ts
    2. If there's a single non-T at the end, check if there are Ts before it
    3. If yes, remove the non-T and continue trimming
    4. Repeat until sequence ends with multiple non-Ts or is all Ts

    Args:
        sequence: The sequence string (uppercase)

    Returns:
        Trimmed sequence
    """
    seq = sequence

    while len(seq) > 0:
        # Remove trailing Ts
        seq_no_t = seq.rstrip('T')

        if len(seq_no_t) == 0:
            # All Ts - return empty or the string based on preference
            return seq_no_t

        if len(seq_no_t) == len(seq):
            # No trailing Ts - we're done
            return seq

        # We removed some Ts. Now check if there's a single non-T at the end
        # If the last character is not T and there are Ts before it, remove it and continue
        if len(seq_no_t) >= 2 and seq_no_t[-1] != 'T' and seq_no_t[-2] == 'T':
            # Pattern like ...TN (where N is the last non-T)
            # Remove the last non-T and continue
            seq = seq_no_t[:-1]
        else:
            # Either ends with multiple non-Ts or just one non-T with no Ts before
            # We're done
            return seq_no_t

    return seq

def run_noeCR34335_pipeline(species_name: str, config: Dict, output_file: Optional[str] = None, strict: bool = False):
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

    # Clean up existing files in results directory for fresh start
    print(f"  Cleaning results directory: {results_dir}", file=sys.stderr)
    for file_path in results_dir.glob("*"):
        if file_path.is_file():
            file_path.unlink()
            print(f"  Removed existing file: {file_path.name}", file=sys.stderr)

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
    extended_file = results_dir / 'lncrna.extended'

    # Get max_polyt_extension from config, default to 2000
    max_extension = config.get('max_polyt_extension', 2000)

    extract_cmd = f"./scripts/extract_with_polyt_extension.py -g {genome_file} -s {species_name} --incomplete {incomplete_file} --extended {extended_file} --max-extension {max_extension} < {filtered_gff} > {temp_output}"

    result = subprocess.run(extract_cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error extracting sequences: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    # Show extension details from stderr
    if result.stderr:
        stderr_lines = result.stderr.split('\n')
        extension_found = False
        i = 0
        while i < len(stderr_lines):
            line = stderr_lines[i].strip()
            if 'Extended' in line:
                if not extension_found:
                    print("  Extension details:", file=sys.stderr)
                    extension_found = True
                print(f"    {line}", file=sys.stderr)
                # Check if next line is "Added sequence"
                if i + 1 < len(stderr_lines) and 'Added sequence' in stderr_lines[i + 1]:
                    added_line = stderr_lines[i + 1].strip()
                    # Truncate long sequences for readability
                    if len(added_line) > 120:
                        sequence_part = added_line.split(': ', 1)[1] if ': ' in added_line else added_line
                        if len(sequence_part) > 80:
                            truncated = sequence_part[:40] + '...' + sequence_part[-40:]
                            print(f"    Added sequence: {truncated} ({len(sequence_part)} bp)", file=sys.stderr)
                        else:
                            print(f"    {added_line}", file=sys.stderr)
                    else:
                        print(f"    {added_line}", file=sys.stderr)
                    i += 1  # Skip the added sequence line since we processed it
            i += 1

    # Use the blast_to_id_mapping created during filtering
    if 'blast_to_id_mapping' not in locals():
        blast_to_id_mapping = {}

    print(f"  Using {len(blast_to_id_mapping)} coordinate mappings for ID replacement", file=sys.stderr)

    # Apply strict filtering if requested
    strict_excluded_ids = set()  # Track IDs excluded by strict filtering
    if strict:
        # Load motifs from config
        motifs_5prime = config.get('motifs', {}).get('5prime', ['GCGGT', 'GTGGT'])
        motifs_3prime = config.get('motifs', {}).get('3prime', ['ATCGC', 'ACCGC', 'ATCAC'])
        motifs_str = '/'.join(motifs_3prime)

        print(f"  Applying strict filtering: excluding sequences without compatible motif pairs and incomplete sequences...", file=sys.stderr)
        # Read extended sequence IDs that need validation
        extended_sequences = set()
        if extended_file.exists():
            with open(extended_file, 'r') as f:
                for line in f:
                    seq_id = line.strip()
                    if seq_id:
                        extended_sequences.add(seq_id)

        # Read incomplete sequence IDs that should be excluded in strict mode
        incomplete_sequences = set()
        if incomplete_file.exists():
            with open(incomplete_file, 'r') as f:
                for line in f:
                    seq_id = line.strip()
                    if seq_id:
                        incomplete_sequences.add(seq_id)

        # Create a temporary file for strict filtering
        temp_strict_output = tempfile.mktemp(suffix='.fa')
        strict_filtered_count = 0
        strict_excluded_count = 0

        with open(temp_output, 'r') as infile, open(temp_strict_output, 'w') as outfile:
            current_header = None
            current_sequence = []

            for line in infile:
                if line.startswith('>'):
                    # Process previous sequence if any
                    if current_header and current_sequence:
                        sequence = ''.join(current_sequence)
                        should_keep = True

                        # Check if this is an extended sequence
                        header_id = current_header[1:].split('|')[0]  # Extract ID from header

                        # Exclude incomplete sequences in strict mode
                        if header_id in incomplete_sequences:
                            should_keep = False
                            strict_excluded_count += 1
                            # Add old ID to exclusion set (will be mapped later)
                            strict_excluded_ids.add(header_id)
                            print(f"    Excluded {header_id}: incomplete sequence (no polyT termination)", file=sys.stderr)

                        # Check motif compatibility for ALL sequences
                        if should_keep:
                            seq_upper = sequence.upper()

                            # Check for 5' motif in first 10 nt
                            first_10nt = seq_upper[:min(10, len(seq_upper))]
                            motif_5prime_found = None
                            for motif_5 in motifs_5prime:
                                if motif_5 in first_10nt:
                                    motif_5prime_found = motif_5
                                    break

                            # For extended sequences, trim polyT properly
                            if header_id in extended_sequences:
                                seq_no_polyt = trim_polyt_proper(seq_upper)
                                seq_for_motif_check = seq_no_polyt + 'T'
                            else:
                                # For non-extended sequences, just trim trailing T's
                                seq_for_motif_check = seq_upper.rstrip('T') + 'T'

                            # Check last 12nt for any of the configured 3' motifs
                            last_12nt = seq_for_motif_check[-12:] if len(seq_for_motif_check) >= 12 else seq_for_motif_check
                            motif_3prime_found = None
                            for motif_3 in motifs_3prime:
                                if motif_3 in last_12nt:
                                    motif_3prime_found = motif_3
                                    break

                            # If 5' motif is found, 3' motif must also be found and compatible
                            if motif_5prime_found:
                                if not motif_3prime_found:
                                    should_keep = False
                                    strict_excluded_count += 1
                                    strict_excluded_ids.add(header_id)
                                    print(f"    Excluded {header_id}: has 5' motif {motif_5prime_found} but no compatible 3' motif in last 12nt: {last_12nt}", file=sys.stderr)
                                else:
                                    # Check motif compatibility
                                    compatibility = config.get('motifs', {}).get('compatibility', {})
                                    compatible_3prime = compatibility.get(motif_5prime_found, [])
                                    if motif_3prime_found not in compatible_3prime:
                                        should_keep = False
                                        strict_excluded_count += 1
                                        strict_excluded_ids.add(header_id)
                                        print(f"    Excluded {header_id}: incompatible motif pair {motif_5prime_found} + {motif_3prime_found}", file=sys.stderr)

                        if should_keep:
                            outfile.write(current_header + '\n')
                            outfile.write(sequence + '\n')
                            strict_filtered_count += 1

                    # Start new sequence
                    current_header = line.strip()
                    current_sequence = []
                else:
                    current_sequence.append(line.strip())

            # Process last sequence
            if current_header and current_sequence:
                sequence = ''.join(current_sequence)
                should_keep = True

                # Check if this is an extended sequence
                header_id = current_header[1:].split('|')[0]  # Extract ID from header

                # Exclude incomplete sequences in strict mode
                if header_id in incomplete_sequences:
                    should_keep = False
                    strict_excluded_count += 1
                    # Add old ID to exclusion set (will be mapped later)
                    strict_excluded_ids.add(header_id)
                    print(f"    Excluded {header_id}: incomplete sequence (no polyT termination)", file=sys.stderr)

                # Check motif compatibility for ALL sequences
                if should_keep:
                    seq_upper = sequence.upper()

                    # Check for 5' motif in first 10 nt
                    first_10nt = seq_upper[:min(10, len(seq_upper))]
                    motif_5prime_found = None
                    for motif_5 in motifs_5prime:
                        if motif_5 in first_10nt:
                            motif_5prime_found = motif_5
                            break

                    # For extended sequences, trim polyT properly
                    if header_id in extended_sequences:
                        seq_no_polyt = trim_polyt_proper(seq_upper)
                        seq_for_motif_check = seq_no_polyt + 'T'
                    else:
                        # For non-extended sequences, just trim trailing T's
                        seq_for_motif_check = seq_upper.rstrip('T') + 'T'

                    # Check last 12nt for any of the configured 3' motifs
                    last_12nt = seq_for_motif_check[-12:] if len(seq_for_motif_check) >= 12 else seq_for_motif_check
                    motif_3prime_found = None
                    for motif_3 in motifs_3prime:
                        if motif_3 in last_12nt:
                            motif_3prime_found = motif_3
                            break

                    # If 5' motif is found, 3' motif must also be found and compatible
                    if motif_5prime_found:
                        if not motif_3prime_found:
                            should_keep = False
                            strict_excluded_count += 1
                            strict_excluded_ids.add(header_id)
                            print(f"    Excluded {header_id}: has 5' motif {motif_5prime_found} but no compatible 3' motif in last 12nt: {last_12nt}", file=sys.stderr)
                        else:
                            # Check motif compatibility
                            compatibility = config.get('motifs', {}).get('compatibility', {})
                            compatible_3prime = compatibility.get(motif_5prime_found, [])
                            if motif_3prime_found not in compatible_3prime:
                                should_keep = False
                                strict_excluded_count += 1
                                strict_excluded_ids.add(header_id)
                                print(f"    Excluded {header_id}: incompatible motif pair {motif_5prime_found} + {motif_3prime_found}", file=sys.stderr)

                if should_keep:
                    outfile.write(current_header + '\n')
                    outfile.write(sequence + '\n')
                    strict_filtered_count += 1

        print(f"  Strict filtering: kept {strict_filtered_count}, excluded {strict_excluded_count} extended sequences", file=sys.stderr)

        # Replace temp_output with strict filtered version
        try:
            os.unlink(temp_output)
        except:
            pass
        temp_output = temp_strict_output

    # Build mapping from old IDs to new IDs while rewriting the file
    old_to_new_id_mapping = {}

    with open(temp_output, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Parse header to extract coordinates
                # Format: >species-number|gene|chrom|start|end|strand
                header = line.strip()
                if '|' in header:
                    parts = header[1:].split('|')
                    if len(parts) >= 6:
                        old_id = parts[0]  # Save old ID
                        chrom = parts[2]
                        start = int(parts[3])
                        end = int(parts[4])
                        strand = parts[5]

                        # Find matching ID by checking for overlapping coordinates
                        # (polyT extension may have changed the coordinates)
                        matched_id = None
                        for original_coord_key, mapped_id in blast_to_id_mapping.items():
                            # Parse original coordinates: "chrom:start-end(strand)"
                            try:
                                coord_parts = original_coord_key.split(':')
                                if len(coord_parts) == 2 and coord_parts[0] == chrom:
                                    range_strand = coord_parts[1].split('(')
                                    if len(range_strand) == 2 and range_strand[1].rstrip(')') == strand:
                                        orig_range = range_strand[0].split('-')
                                        if len(orig_range) == 2:
                                            orig_start = int(orig_range[0])
                                            orig_end = int(orig_range[1])

                                            # Check if coordinates overlap or are close (within 1000bp extension range)
                                            # For - strand, extension goes upstream (start decreases)
                                            # For + strand, extension goes downstream (end increases)
                                            overlap = False
                                            if strand == '-':
                                                # Check if extended start <= original start and end matches or overlaps
                                                if start <= orig_start and abs(end - orig_end) <= 50:
                                                    overlap = True
                                            else:
                                                # Check if start matches or overlaps and extended end >= original end
                                                if abs(start - orig_start) <= 50 and end >= orig_end:
                                                    overlap = True

                                            if overlap:
                                                matched_id = mapped_id
                                                break
                            except (ValueError, IndexError):
                                continue

                        if matched_id:
                            old_to_new_id_mapping[old_id] = matched_id  # Save mapping
                            outfile.write(f">{matched_id}|noeCR34335|{chrom}|{start}|{end}|{strand}\n")
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
                old_id = line.strip()
                if old_id:
                    # Look up the new ID from the mapping we just created
                    new_id = old_to_new_id_mapping.get(old_id, old_id)
                    updated_incomplete.append(new_id)

        # Write updated incomplete file
        with open(incomplete_file, 'w') as f:
            for seq_id in updated_incomplete:
                f.write(f"{seq_id}\n")
        print(f"  Updated {len(updated_incomplete)} incomplete sequence IDs", file=sys.stderr)

    # Update extended file with mapped IDs if it exists
    if extended_file.exists():
        print(f"  Updating extended sequence IDs...", file=sys.stderr)
        updated_extended = []

        with open(extended_file, 'r') as f:
            for line in f:
                old_id = line.strip()
                if old_id:
                    # Look up the new ID from the mapping we created earlier
                    new_id = old_to_new_id_mapping.get(old_id, old_id)
                    updated_extended.append(new_id)

        # Write updated extended file
        with open(extended_file, 'w') as f:
            for seq_id in updated_extended:
                f.write(f"{seq_id}\n")
        print(f"  Updated {len(updated_extended)} extended sequence IDs", file=sys.stderr)

    # Step 5: Copy filtered upstream sequences to results directory with updated IDs
    upstream_output = results_dir / "upstream.fa"
    print(f"  Copying filtered upstream sequences to {upstream_output}...", file=sys.stderr)

    # Create mapping from coordinates to new IDs
    coord_to_id = {}
    for (chrom, start, end, strand, counter) in hits_with_promoters:
        coord_to_id[(chrom, start, end, strand)] = f"{species_abbrev}{counter}"

    # Write sequences in 2-line format to match original
    upstream_sequences_written = 0
    with open(upstream_output, 'w') as f:
        for record in upstream_with_promoters:
            should_write = True
            header_to_write = None

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
                        # Check if this sequence was excluded by strict filtering
                        # We need to check if any excluded old ID would map to this new_id by checking coordinates
                        is_excluded = False
                        if strict:
                            # Find the ncRNA coordinates that correspond to this upstream region
                            # Upstream is CM025933.1:21374518-21374617(+)
                            # ncRNA is CM025933.1:21374618-21374678(+) (starts right after upstream)
                            for coord_key in blast_to_id_mapping:
                                if blast_to_id_mapping[coord_key] == new_id:
                                    # This ncRNA maps to new_id, now check if any excluded old_id has these coordinates
                                    # We stored old IDs in strict_excluded_ids, but they don't have coordinates
                                    # We need to check the incomplete file which was read before updating IDs
                                    # Actually, simpler: just check if new_id maps to any excluded sequence
                                    # by checking if old_to_new_id_mapping has the reverse mapping
                                    for old_id in strict_excluded_ids:
                                        if old_to_new_id_mapping.get(old_id) == new_id:
                                            is_excluded = True
                                            break
                                    # If old_id not in mapping (because excluded before mapping), check by matching coord
                                    if not is_excluded:
                                        # Check if this coord is for an excluded sequence by checking if
                                        # the sequence with these coords was excluded
                                        # Since excluded sequences aren't in old_to_new_id_mapping, we need another way
                                        # The excluded ID is like "Drosophila_subpulchrella-4" which doesn't have coords
                                        # But we know new_id from coord_to_id mapping
                                        # So if there are N total sequences and M kept, check if new_id counter is in excluded range
                                        # Actually simpler: just check all values in old_to_new_id_mapping
                                        # If new_id is NOT in the values, it means it was excluded
                                        if new_id not in old_to_new_id_mapping.values():
                                            is_excluded = True
                                    break

                        if is_excluded:
                            print(f"    Excluding upstream sequence for {new_id} (sequence was filtered out)", file=sys.stderr)
                            should_write = False
                        else:
                            # Create new header with updated ID but keeping noeCR34335 as gene type
                            header_to_write = f">{new_id}|noeCR34335|{chrom}|{start}|{end}|{strand}\n"
                    else:
                        header_to_write = f">{record.id}\n"
                else:
                    header_to_write = f">{record.id}\n"
            else:
                header_to_write = f">{record.id}\n"

            if should_write:
                f.write(header_to_write)
                f.write(f"{str(record.seq)}\n")
                upstream_sequences_written += 1

    print(f"  Wrote {upstream_sequences_written} upstream sequences with promoters to {upstream_output}", file=sys.stderr)

    # Score noeCR34335 promoters using new PWM scoring script
    score_output = results_dir / "upstream.score"
    print(f"  Scoring noeCR34335 promoters...", file=sys.stderr)

    # Use the new PWM scoring pipeline command
    pwm_result = subprocess.run(f"cat {upstream_output} | ./scripts/score_sequences_with_pwm.py -e noeCR34335 -c config.yaml",
                               shell=True, capture_output=True, text=True)

    if pwm_result.returncode != 0:
        print(f"  Warning: PWM scoring failed: {pwm_result.stderr}", file=sys.stderr)
        result = pwm_result
    else:
        # Extract noeCR34335 lines and filter by Pol III score threshold from config
        filtered_ids = set()  # Track IDs that pass the filter
        pol3_threshold = config.get('pse_score_cutoff', 0.8)

        with open(score_output, 'w') as f:
            f.write("ID\tSequence\tPol2_Score\tPol3_Score\n")

            for line in pwm_result.stdout.strip().split('\n'):
                if 'noeCR34335' in line and '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        full_seq_id = parts[0]
                        pol2_score = parts[1]
                        pol3_score = parts[2]

                        # Check if Pol III score meets threshold
                        try:
                            pol3_value = float(pol3_score)
                            if pol3_value < pol3_threshold:
                                print(f"    Filtering out {full_seq_id}: Pol III score {pol3_score} < {pol3_threshold}", file=sys.stderr)
                                continue  # Skip this sequence
                        except ValueError:
                            print(f"    Warning: Invalid Pol III score '{pol3_score}' for {full_seq_id}", file=sys.stderr)
                            continue

                        # Extract just the short ID for the output
                        seq_id = full_seq_id.split('|')[0]
                        filtered_ids.add(seq_id)

                        # Extract uppercase sequence from original file using full ID
                        seq_cmd = f"grep -A1 '^>{full_seq_id}$' {upstream_output} | tail -1 | tr -d 'a-z-'"
                        seq_result = subprocess.run(seq_cmd, shell=True, capture_output=True, text=True)
                        sequence = seq_result.stdout.strip() if seq_result.returncode == 0 else ""

                        f.write(f"{seq_id}\t{sequence}\t{pol2_score}\t{pol3_score}\n")

        # Filter lncrna.fa and upstream.fa to remove noeCR34335 sequences with low Pol III scores
        print(f"  Filtering FASTA files to remove noeCR34335 sequences with Pol III score < {pol3_threshold}", file=sys.stderr)

        # Filter lncrna.fa (only remove low-scoring noeCR34335 sequences)
        from Bio import SeqIO
        temp_lncrna = output_file.with_suffix('.temp.fa')
        with open(temp_lncrna, 'w') as out_f:
            for record in SeqIO.parse(output_file, 'fasta'):
                seq_id = record.id.split('|')[0]
                # Keep sequence if it's noeCR34335 and passed filter, or if it's not noeCR34335
                if 'noeCR34335' in record.id:
                    if seq_id in filtered_ids:  # Only keep noeCR34335 that passed filter
                        SeqIO.write(record, out_f, 'fasta')
                else:
                    # This shouldn't happen in lncrna.fa, but keep for safety
                    SeqIO.write(record, out_f, 'fasta')

        # Replace original with filtered version
        os.rename(temp_lncrna, output_file)

        # Filter upstream.fa (keep all non-noeCR34335, only filter noeCR34335 by score)
        temp_upstream = upstream_output.with_suffix('.temp.fa')
        with open(temp_upstream, 'w') as out_f:
            for record in SeqIO.parse(upstream_output, 'fasta'):
                seq_id = record.id.split('|')[0]
                # Keep all non-noeCR34335 sequences, and only noeCR34335 that passed filter
                if 'noeCR34335' in record.id:
                    if seq_id in filtered_ids:  # Only keep noeCR34335 that passed filter
                        SeqIO.write(record, out_f, 'fasta')
                else:
                    # Keep all other ncRNA types (U1, U2, U6, etc.)
                    SeqIO.write(record, out_f, 'fasta')

        # Replace original with filtered version
        os.rename(temp_upstream, upstream_output)

        print(f"  Kept {len(filtered_ids)} noeCR34335 sequences and all other ncRNA types", file=sys.stderr)

        result = type('obj', (object,), {'returncode': 0, 'stderr': ''})()

    if result.returncode != 0:
        print(f"  Warning: PWM scoring failed: {result.stderr}", file=sys.stderr)
    else:
        print(f"  Wrote promoter scores to {score_output}", file=sys.stderr)

    # Step 9: Additional filtering - Check if genomic RNA sequences start with G* followed by 5' motif pattern
    # Load 5' motifs from config
    motifs_5prime = config.get('motifs', {}).get('5prime', ['GCGGT', 'GTGGT'])
    motifs_pattern = '|'.join(motifs_5prime)
    print(f"\n  Checking genomic RNA sequences for G*({motifs_pattern}) pattern...", file=sys.stderr)
    import re
    from Bio import SeqIO
    pattern = re.compile(rf'^G*({motifs_pattern})', re.IGNORECASE)

    sequences_to_keep = []
    sequences_removed = []

    # Read and check each sequence in lncrna.fa
    if output_file.exists():
        for record in SeqIO.parse(output_file, 'fasta'):
            seq_str = str(record.seq).upper()
            if pattern.match(seq_str):
                sequences_to_keep.append(record)
            else:
                sequences_removed.append(record.id)
                # Print red warning
                print(f"    \033[91mWARNING: Removing {record.id} - sequence does not match /^G*({motifs_pattern})/ pattern\033[0m", file=sys.stderr)
                print(f"    \033[91m         Sequence starts with: {seq_str[:50]}...\033[0m", file=sys.stderr)

        # Rewrite lncrna.fa with only sequences that match the pattern
        if sequences_removed:
            with open(output_file, 'w') as f:
                SeqIO.write(sequences_to_keep, f, 'fasta')

            # Get short IDs of removed sequences
            removed_short_ids = set()
            for record_id in sequences_removed:
                short_id = record_id.split('|')[0]
                removed_short_ids.add(short_id)

            # Also update upstream.fa to remove corresponding sequences
            if upstream_output.exists():
                upstream_to_keep = []
                for record in SeqIO.parse(upstream_output, 'fasta'):
                    seq_id = record.id.split('|')[0]
                    # Only remove noeCR34335 sequences that failed the pattern test
                    if 'noeCR34335' in record.id and seq_id in removed_short_ids:
                        continue  # Skip this sequence
                    upstream_to_keep.append(record)

                # Rewrite upstream.fa
                with open(upstream_output, 'w') as f:
                    SeqIO.write(upstream_to_keep, f, 'fasta')

            print(f"  Filtered out {len(sequences_removed)} sequences not matching G*({motifs_pattern}) pattern", file=sys.stderr)
            print(f"  Kept {len(sequences_to_keep)} sequences in lncrna.fa", file=sys.stderr)
        else:
            print(f"  All sequences match the G*({motifs_pattern}) pattern", file=sys.stderr)

    # Clean up temporary files
    try:
        os.unlink(temp_gff)
        os.unlink(filtered_gff)
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

    # Step 8: Calculate pairwise identities if there are multiple sequences
    identities_file = results_dir / "lncrna.identities"
    if output_file.exists():
        # Count sequences in lncrna.fa
        seq_count = 0
        with open(output_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    seq_count += 1

        if seq_count >= 2:
            print(f"\n  Calculating pairwise identities for {seq_count} sequences...", file=sys.stderr)
            identity_cmd = f"cat {output_file} | ./scripts/calculate_pairwise_identities.py > {identities_file}"
            result = subprocess.run(identity_cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"  Warning: Pairwise identity calculation failed: {result.stderr}", file=sys.stderr)
                # Create empty file to avoid errors later
                identities_file.touch()
            else:
                print(f"  Pairwise identities written to {identities_file}", file=sys.stderr)
        else:
            print(f"  Skipping pairwise identity calculation (only {seq_count} sequence(s))", file=sys.stderr)
            # Create empty file
            identities_file.touch()
    else:
        # Create empty file if lncrna.fa doesn't exist
        identities_file.touch()

    # Step 9: Calculate sequence statistics
    stats_file = results_dir / "lncrna.stats"
    if output_file.exists():
        print(f"\n  Calculating sequence statistics...", file=sys.stderr)
        stats_cmd = f"./scripts/calculate_lncrna_stats.py {output_file} -o {stats_file}"
        result = subprocess.run(stats_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  Warning: Statistics calculation failed: {result.stderr}", file=sys.stderr)
            # Create empty file to avoid errors later
            stats_file.touch()
        else:
            print(f"  Sequence statistics written to {stats_file}", file=sys.stderr)
    else:
        # Create empty file if lncrna.fa doesn't exist
        stats_file.touch()

    # Step 10: Generate HTML and PDF reports
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
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>lncRNA Analysis</title>
    <style>
        @page {
            size: A4 portrait;
            margin: 1.5cm;
        }

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: 'Arial', 'Helvetica', sans-serif;
            font-size: 9pt;
            line-height: 1.3;
            color: #000;
            background: white;
            max-width: 210mm;
            margin: 0 auto;
            padding: 10px;
        }

        /* Compact headers */
        h1 {
            font-size: 12pt;
            font-weight: bold;
            margin-bottom: 4px;
            font-style: italic;
            border-bottom: 1px solid #333;
            padding-bottom: 2px;
        }

        h2 {
            font-size: 10pt;
            font-weight: bold;
            margin-top: 8px;
            margin-bottom: 4px;
            color: #222;
        }

        /* Assembly info - inline */
        .assembly-info {
            font-size: 8pt;
            margin-bottom: 8px;
            color: #444;
        }

        /* Taxonomic lineage */
        .taxonomic-lineage {
            font-size: 8pt;
            margin-bottom: 6px;
            color: #444;
        }

        /* Compact sequence container */
        .sequence-entry {
            border: 1px solid #ccc;
            padding: 4px 6px;
            margin: 6px 0;
            background: #fafafa;
        }

        /* Sequence header - single line */
        .sequence-header {
            font-size: 8pt;
            font-weight: bold;
            margin-bottom: 2px;
        }

        /* Inline metadata */
        .sequence-metadata {
            font-size: 7pt;
            color: #666;
            margin-bottom: 3px;
            line-height: 1.4;
        }

        .metadata-item {
            display: inline-block;
            margin-right: 8px;
            white-space: nowrap;
        }

        /* Compact score labels */
        .score {
            display: inline;
            padding: 0 2px;
            font-weight: bold;
            border-radius: 2px;
        }

        .score-high {
            background-color: #28a745;
            color: white;
        }

        .score-medium {
            background-color: #ffc107;
            color: #333;
        }

        .score-low {
            background-color: #dc3545;
            color: white;
        }

        .status-ok {
            color: #28a745;
            font-weight: bold;
        }

        /* Sequence display */
        .sequence {
            font-family: 'Courier New', Courier, monospace;
            font-size: 5pt;
            letter-spacing: 0.02em;
            line-height: 1.0;
            margin: 2px 0;
        }

        /* Compact alignment */
        .alignment-section {
            margin-top: 10px;
            border-top: 1px solid #ccc;
            padding-top: 6px;
        }

        .alignment {
            font-family: 'Courier New', Courier, monospace;
            font-size: 5pt;
            line-height: 1.0;
            letter-spacing: 0.01em;
            background: #fafafa;
            border: 1px solid #ddd;
            padding: 4px 6px;
            margin: 4px 0;
            overflow-x: auto;
            white-space: pre-wrap;
        }

        .alignment-row {
            margin-bottom: 1px;
        }

        .alignment-label {
            display: inline-block;
            width: 90px;
            font-weight: bold;
            color: #333;
            font-size: 7pt;
        }

        /* Separator between entries */
        .entry-separator {
            border-bottom: 1px dotted #999;
            margin: 8px 0;
        }

        /* Print optimization */
        @media print {
            body {
                font-size: 8pt;
                padding: 0;
            }
            .sequence-entry {
                page-break-inside: avoid;
            }
            .alignment {
                page-break-inside: avoid;
            }
        }

        /* Bold sequences styling - removed underline */
    </style>
</head>
<body>
""")

        # Species information
        species_display = species_name.replace('_', ' ')
        assembly_name = get_assembly_name(species_name)
        f.write(f"<h1>{species_display}</h1>\n")

        # Get taxonomic lineage
        lineage_cmd = f"python3 ./scripts/get_taxonomic_lineage.py '{species_display}'"
        lineage_result = subprocess.run(lineage_cmd, shell=True, capture_output=True, text=True)
        if lineage_result.returncode == 0 and lineage_result.stdout.strip() and not lineage_result.stdout.startswith("Taxonomic lineage not found"):
            # Remove the species name from the end of the lineage
            lineage = lineage_result.stdout.strip()
            lineage_parts = lineage.split(' > ')
            # Remove the last part if it matches the species name
            if lineage_parts and lineage_parts[-1].lower() == species_display.lower().replace('_', ' '):
                lineage_parts = lineage_parts[:-1]
            lineage_without_species = ' > '.join(lineage_parts)

            f.write(f'<div class="taxonomic-lineage">\n')
            f.write(f"<strong>Taxonomic lineage:</strong> {lineage_without_species}\n")
            f.write("</div>\n")

        f.write(f'<div class="assembly-info">\n')
        f.write(f"<strong>Assembly:</strong> {assembly_name}\n")
        f.write("</div>\n")

        # lncRNA sequences section
        # Format lncRNA sequences with scores
        if output_file.exists():
            format_cmd = f"cat {output_file} | python3 ./scripts/html_format_seq.py {score_output} {incomplete_file} {noe_scores_file} {stats_file} {extended_file}"
            result = subprocess.run(format_cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                f.write(result.stdout)
            else:
                f.write("Error formatting lncRNA sequences\n")

        # Pairwise sequence identities section - AFTER sequences, BEFORE upstream alignment
        if identities_file.exists() and identities_file.stat().st_size > 0:
            # Read and parse the identities file
            identities_data = []
            seq_ids = []
            identity_dict = {}

            with open(identities_file, 'r') as id_f:
                header = id_f.readline()  # Skip header
                for line in id_f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            seq1, seq2, identity_val = parts[0], parts[1], float(parts[2])
                            identities_data.append((seq1, seq2, identity_val))

                            # Track unique sequence IDs in order
                            if seq1 not in seq_ids:
                                seq_ids.append(seq1)
                            if seq2 not in seq_ids:
                                seq_ids.append(seq2)

                            # Store in dictionary for matrix lookup
                            identity_dict[(seq1, seq2)] = identity_val
                            identity_dict[(seq2, seq1)] = identity_val

            if identities_data and len(seq_ids) <= 10:
                # For 10 or fewer sequences, show a compact matrix
                f.write('<div style="margin: 10px 0;">\n')
                f.write('<table style="border-collapse:collapse;font-size:9px;">\n')

                # Header row
                f.write('<tr><td style="border:1px solid #ddd;padding:3px;background:#f2f2f2;"></td>')
                for seq_id in seq_ids:
                    f.write(f'<th style="border:1px solid #ddd;padding:3px;background:#f2f2f2;font-weight:bold;">{seq_id}</th>')
                f.write('</tr>\n')

                # Data rows
                for i, seq1 in enumerate(seq_ids):
                    f.write(f'<tr><th style="border:1px solid #ddd;padding:3px;background:#f2f2f2;text-align:left;font-weight:bold;">{seq1}</th>')
                    for j, seq2 in enumerate(seq_ids):
                        if i == j:
                            # Diagonal - self comparison
                            f.write('<td style="border:1px solid #ddd;padding:3px;text-align:center;background:#e0e0e0;">-</td>')
                        elif i < j:
                            # Upper triangle - show values
                            identity = identity_dict.get((seq1, seq2), None)
                            if identity is None:
                                # Missing comparison
                                f.write('<td style="border:1px solid #ddd;padding:3px;text-align:center;background:#f0f0f0;">n/a</td>')
                            else:
                                # Show identity value without color coding
                                f.write(f'<td style="border:1px solid #ddd;padding:3px;text-align:center;font-size:8px;">')
                                f.write(f'{identity:.0f}</td>')
                        else:
                            # Lower triangle - empty
                            f.write('<td style="border:1px solid #ddd;padding:3px;background:#f9f9f9;"></td>')
                    f.write('</tr>\n')

                f.write('</table>\n')
                f.write('<div style="font-size:7pt;margin-top:3px;color:#333;">\n')
                f.write('Pairwise sequence identity matrix (%) showing similarity between lncRNA sequences.\n')
                f.write('</div>\n')
                f.write('</div>\n')

            elif identities_data:
                # For more than 10 sequences, show summary statistics only
                all_identities = [identity for _, _, identity in identities_data]
                min_identity = min(all_identities)
                max_identity = max(all_identities)
                avg_identity = sum(all_identities) / len(all_identities)

                # Calculate distribution
                low_count = sum(1 for i in all_identities if i < 30)
                med_count = sum(1 for i in all_identities if 30 <= i < 60)
                high_count = sum(1 for i in all_identities if i >= 60)

                f.write('<div style="margin: 10px 0;">\n')
                f.write('<h3 style="font-size:11pt;margin-bottom:5px;">Pairwise Identity Summary</h3>\n')
                f.write('<div style="background:#f9f9f9;padding:8px;border:1px solid #ddd;border-radius:3px;font-size:10px;">\n')
                f.write(f'<div><strong>Sequences:</strong> {len(seq_ids)}</div>\n')
                f.write(f'<div><strong>Comparisons:</strong> {len(identities_data)}</div>\n')
                f.write(f'<div><strong>Identity range:</strong> {min_identity:.1f}% - {max_identity:.1f}% (avg: {avg_identity:.1f}%)</div>\n')
                f.write('<div><strong>Distribution:</strong> ')
                f.write(f'Low (&lt;30%): {low_count} | ')
                f.write(f'Medium (30-60%): {med_count} | ')
                f.write(f'High (≥60%): {high_count}')
                f.write('</div>\n')
                f.write('</div>\n')
                f.write('</div>\n')

        # Upstream promoter alignment section
        f.write('<div class="alignment-section">\n')
        f.write('<pre class="alignment">\n')

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

        # Generate description of Pol II and Pol III types from config
        pol2_types = config.get('rna_polymerases', {}).get('pol2', {}).get('types', [])
        pol3_types = config.get('rna_polymerases', {}).get('pol3', {}).get('types', [])

        f.write('<div style="font-size:7pt;margin-top:3px;color:#333;">\n')
        f.write('Upstream promoter alignment showing conserved regions. ')
        if pol3_types:
            f.write(f'Sequences include Pol III promoters ({", ".join(pol3_types)}) ')
        if pol2_types:
            f.write(f'and Pol II promoters ({", ".join(pol2_types)}). ')
        f.write('Conserved positions are highlighted in the alignment.\n')
        f.write('</div>\n')
        f.write("</div>\n")

        f.write("</body>\n</html>\n")

    print(f"  Generated HTML report: {html_file}", file=sys.stderr)

    # Generate PDF using headless Chrome
    pdf_cmd = f"google-chrome --headless --disable-gpu --no-pdf-header-footer --print-to-pdf={pdf_file} {html_file}"
    result = subprocess.run(pdf_cmd, shell=True, capture_output=True, text=True)

    if result.returncode == 0 and pdf_file.exists():
        print(f"  Generated PDF report: {pdf_file}", file=sys.stderr)
    else:
        print(f"  Warning: PDF generation failed. HTML report available at {html_file}", file=sys.stderr)

    # Check if lncrna.fa is empty and clean up if needed
    if output_file.exists():
        # Count sequences in the FASTA file
        sequence_count = 0
        try:
            for record in SeqIO.parse(output_file, 'fasta'):
                sequence_count += 1
                break  # We only need to know if there's at least one sequence
        except:
            sequence_count = 0

        if sequence_count == 0:
            print(f"  Warning: No lncRNA sequences found after filtering", file=sys.stderr)
            print(f"  Removing empty results directory: {results_dir}", file=sys.stderr)
            import shutil
            shutil.rmtree(results_dir)
            print(f"noeCR34335 lncRNA extraction completed for {species_name} - NO RESULTS (directory removed)", file=sys.stderr)
        else:
            print(f"\nnoeCR34335 lncRNA extraction completed for {species_name}", file=sys.stderr)
    else:
        print(f"  Error: Output file {output_file} was not created", file=sys.stderr)
        print(f"  Removing results directory: {results_dir}", file=sys.stderr)
        import shutil
        shutil.rmtree(results_dir)
        print(f"noeCR34335 lncRNA extraction completed for {species_name} - NO RESULTS (directory removed)", file=sys.stderr)

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
    parser.add_argument('--strict', action='store_true',
                       help='Filter out extended sequences without ATCGC motif (for noeCR34335 command)')

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
                output_file=args.output,
                strict=args.strict
            )

    if len(commands) > 1:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"PIPELINE COMPLETE: All {len(commands)} steps finished successfully!", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)

if __name__ == "__main__":
    main()