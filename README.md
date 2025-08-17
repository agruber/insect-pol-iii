# Insect Pol-III type 3

Bioinformatics research project focused on FASTA sequence manipulation and analysis using command-line tools and BioPython.

## Genome Download and Search Pipeline

Unix-style pipeline for bioinformatics workflows with three core scripts:
1. **filter_species.py** - Filter species based on taxonomic queries
2. **download_genomes.py** - Generate download commands (reads from stdin)
3. **search_genomes.py** - Generate cmsearch commands (reads from stdin)

### Pipeline Architecture

The scripts follow Unix philosophy - each does one thing well and can be chained together:

```
filter_species.py → download_genomes.py
filter_species.py → search_genomes.py
```

### Scripts

#### 1. Species Filter (`filter_species.py`)

Core filtering script that outputs species data based on taxonomic queries.

**Usage:**
```bash
# Species search (TSV output)
./scripts/filter_species.py data/genomes_annotated.tsv "Drosophila melanogaster"

# Species search (names only)
./scripts/filter_species.py data/genomes_annotated.tsv "Drosophila melanogaster" -f names

# Genus search
./scripts/filter_species.py data/genomes_annotated.tsv "Drosophila" -t genus

# Order search
./scripts/filter_species.py data/genomes_annotated.tsv "Lepidoptera" -t order
```

**Output Formats:**
- `tsv`: Full TSV data (default) - for download pipeline
- `names`: Cleaned species names only - for search pipeline

#### 2. Genome Downloader (`download_genomes.py`)

Reads TSV data from stdin and generates download commands.

**Usage:**
```bash
# Download single species
./scripts/filter_species.py data/genomes_annotated.tsv "Drosophila melanogaster" | ./scripts/download_genomes.py

# Download genus
./scripts/filter_species.py data/genomes_annotated.tsv "Drosophila" -t genus | ./scripts/download_genomes.py
```

**Output Structure:**
```
genomes/
├── Species_name/
│   ├── genome.fna.gz -> actual_file.fna.gz  # Symlink
│   ├── actual_file.fna.gz                   # Downloaded genome
│   └── logs/
│       └── download.txt                     # Download metadata
```

#### 3. cmsearch Runner (`search_genomes.py`)

Reads species names from stdin and generates cmsearch commands.

**Usage:**
```bash
# Search single species
./scripts/filter_species.py data/genomes_annotated.tsv "Drosophila melanogaster" -f names | ./scripts/search_genomes.py data/models/all.cm

# Search with batching
./scripts/filter_species.py data/genomes_annotated.tsv "Drosophila" -t genus -f names | ./scripts/search_genomes.py data/models/all.cm -b 10
```

**Features:**
- Checks if genome exists before processing
- Skips already processed species
- Decompresses genome.fna.gz to temp file
- Runs cmsearch with --tblout and --noali
- Compresses results and cleans up temp files
- Supports batching for large datasets

**Output:**
- Results stored as `cmsearch_results.txt.gz` in each species directory

### Workflow Examples

#### Example: Drosophila melanogaster

1. **Download genome:**
```bash
./scripts/filter_species.py data/genomes_annotated.tsv "Drosophila melanogaster" | ./scripts/download_genomes.py | bash
```

2. **Run cmsearch:**
```bash
./scripts/filter_species.py data/genomes_annotated.tsv "Drosophila melanogaster" -f names | ./scripts/search_genomes.py data/models/all.cm | bash
```

#### Large-scale Workflow

1. **Download genomes for a taxonomic group:**
```bash
./scripts/filter_species.py data/genomes_annotated.tsv "Lepidoptera" -t order | ./scripts/download_genomes.py > download_lepidoptera.sh
chmod +x download_lepidoptera.sh
./download_lepidoptera.sh
```

2. **Run cmsearch with batching:**
```bash
./scripts/filter_species.py data/genomes_annotated.tsv "Lepidoptera" -t order -f names | ./scripts/search_genomes.py data/models/all.cm -b 20 > search_lepidoptera.sh
chmod +x search_lepidoptera.sh
./search_lepidoptera.sh
```

#### Batch Processing

For large datasets, split work into manageable batches:

```bash
# Generate search commands with batches of 10
./scripts/filter_species.py data/genomes_annotated.tsv "Insecta" -t class -f names | ./scripts/search_genomes.py data/models/all.cm -b 10 > search_insecta.sh

# Split into separate batch files
csplit search_insecta.sh '/# Batch [0-9]/' '{*}' --prefix=batch_ --suffix-format=%02d.sh

# Run batches in parallel
for batch in batch_*.sh; do
    chmod +x "$batch"
    nohup "./$batch" > "${batch%.sh}.log" 2>&1 &
done
```

### Taxonomic Levels Supported

- phylum, class, order, suborder, infraorder
- superfamily, family, genus, subgenus
- species_group, species_subgroup

### Requirements

- Python 3 with pandas
- wget (for downloads)
- cmsearch (INFERNAL)
- gunzip/gzip
- Models file at `data/models/all.cm`
- Genomes metadata at `data/genomes_annotated.tsv`

## Data Acquisition and Processing

### 1. Download Open Tree of Life Taxonomy

Download the taxonomy database from the Open Tree of Life project:

```bash
wget -O data/ott3.7.2.tgz https://files.opentreeoflife.org/ott/ott3.7/ott3.7.2.tgz 
tar -xzf data/ott3.7.2.tgz -C data/ --strip-components=1 ott3.7.2/taxonomy.tsv
```

This provides comprehensive taxonomic information for species annotation.

### 2. Download Invertebrate Genome Data

Retrieve representative and reference genomes from NCBI GenBank:

```bash
curl https://ftp.ncbi.nih.gov/genomes/genbank/invertebrate/assembly_summary.txt | \
  grep -P '(representative genome|reference genome)' | \
  cut -f8,12,20 | sort > data/genomes.tsv
```

This creates a TSV file with:
- Column 1: Species name
- Column 2: Assembly level (Chromosome/Scaffold/Contig)
- Column 3: FTP download URL

### 3. Add Taxonomic Annotations

Use the optimized taxonomy annotation script to add taxonomic information:

```bash
./scripts/add_taxonomy.py data/genomes.tsv data/taxonomy.tsv -v -o results/annotated_genomes.tsv
```

**Output includes 11 additional taxonomic columns:**
- Phylum
- Class  
- Order
- Suborder
- Infraorder
- Superfamily
- Family
- Genus
- Subgenus
- Species group
- Species subgroup
