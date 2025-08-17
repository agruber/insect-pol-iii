# Insect Pol-III type 3

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
