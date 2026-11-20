#!/bin/bash
# Downloading RNA fastq Files from ENA & ATAC Files Pre-fetching from PRJNA998429

cd /COMMA/data/meta_analysis/fastq/multiome_leblanc

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/002/SRR25421402/SRR25421402_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/003/SRR25421403/SRR25421403_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/005/SRR25421405/SRR25421405_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/000/SRR25421400/SRR25421400_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/097/SRR25421397/SRR25421397_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/095/SRR25421395/SRR25421395_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/004/SRR25421404/SRR25421404_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/094/SRR25421394/SRR25421394_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/099/SRR25421399/SRR25421399_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/092/SRR25421392/SRR25421392_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/095/SRR25421395/SRR25421395_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/096/SRR25421396/SRR25421396_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/098/SRR25421398/SRR25421398_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/099/SRR25421399/SRR25421399_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/092/SRR25421392/SRR25421392_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/001/SRR25421401/SRR25421401_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/005/SRR25421405/SRR25421405_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/003/SRR25421403/SRR25421403_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/093/SRR25421393/SRR25421393_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/093/SRR25421393/SRR25421393_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/001/SRR25421401/SRR25421401_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/096/SRR25421396/SRR25421396_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/002/SRR25421402/SRR25421402_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/097/SRR25421397/SRR25421397_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/098/SRR25421398/SRR25421398_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/004/SRR25421404/SRR25421404_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/094/SRR25421394/SRR25421394_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/000/SRR25421400/SRR25421400_1.fastq.gz

export PATH="$PATH:software/sratoolkit.3.1.1-centos_linux64/bin"
while read -r file; do
  prefetch $file --max-size 100g 
done < /COMMA/data/meta_analysis/fastq/multiome_leblanc/accession_list.txt