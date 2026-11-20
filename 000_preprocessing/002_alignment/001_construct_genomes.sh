#!/bin/bash

export PATH="$PATH:software/cellranger-8.0.1/bin"

module load STAR/2.7.11b-GCC-12.3.0

# Downloading Primary Assembly & Basic Annotations from GENCODE

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz 

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz
gunzip gencode.v46.basic.annotation.gtf.gz

cd /COMMA/reference_genome/GENCODE/GRCh38.p14.GENCODE_release46

# Filtering GTF Files for Protein Coding and lncRNA Genes

cellranger mkgtf gencode.v46.basic.annotation.gtf \
  gencode.v46.basic.annotation_filtered.gtf \
  --attribute=gene_type:protein_coding \
  --attribute=gene_type:lncRNA 

cd /COMMA/reference_genome/

# Consturcting STAR Genome

STAR  --runMode genomeGenerate --runThreadN 6 --genomeDir ./STAR_meta --genomeFastaFiles '/COMMA/reference_genome/GENCODE/GRCh38.p14.GENCODE_release46/GRCh38.primary_assembly.genome.fa'  --sjdbGTFfile '/COMMA/reference_genome/GENCODE/GRCh38.p14.GENCODE_release46/gencode.v46.basic.annotation_filtered.gtf'

# Constructing CellRanger Genome

cellranger mkref \
  --genome GRCh38 \
  --fasta '/COMMA/reference_genome/GENCODE/GRCh38.p14.GENCODE_release46/GRCh38.primary_assembly.genome.fa' \
  --genes '/COMMA/reference_genome/GENCODE/GRCh38.p14.GENCODE_release46/gencode.v46.basic.annotation_filtered.gtf'

# Downloading TF Motif File
wget https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt
export PATH="$PATH:software/cellranger-arc-2.0.2/bin"

# Constructing Genome for CellRanger ARC and ATAC

cd /COMMA/reference_genome/cellranger_meta/cellranger-arc_gencode46_jaspar2024
cellranger-arc mkref --config /COMMA/reference_genome/cellranger_meta/cellranger-arc_gencode46_jaspar2024/config