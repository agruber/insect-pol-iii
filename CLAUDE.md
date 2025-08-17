# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview
Bioinformatics research project focused on FASTA sequence manipulation and analysis using command-line tools and BioPython.

## Environment Setup
- **Python**: 3.8.10 with BioPython 1.81, pandas, numpy, matplotlib, scipy
- **Missing tools**: seqtk, samtools, bedtools, BLAST+, fastp
- **Package manager**: conda/mamba recommended for bioinformatics tools

## Installation Commands
```bash
# Install essential bioinformatics tools
conda install -c bioconda seqtk samtools bedtools blast fastp

# Or using apt
sudo apt install seqtk samtools bedtools ncbi-blast+ fastp

# Additional Python packages
pip3 install seaborn jupyter plotly
```

## Project Structure
```
├── scripts/          # Analysis scripts and utilities
├── data/
│   ├── raw/         # Original input data
│   └── processed/   # Cleaned/processed data
├── results/         # Analysis outputs
└── notebooks/       # Jupyter notebooks for interactive analysis
```

## Key Tools for FASTA Manipulation
- **seqtk**: Swiss Army knife for FASTA/FASTQ operations
- **samtools faidx**: FASTA indexing and extraction
- **bedtools**: Genomic interval operations
- **BioPython**: Programmatic sequence manipulation
- **BLAST+**: Sequence similarity searches

## Common Workflows
- Use command-line tools for bulk operations and preprocessing
- Use BioPython for complex programmatic analysis
- Store intermediate results in `data/processed/`
- Use Jupyter notebooks for exploratory analysis