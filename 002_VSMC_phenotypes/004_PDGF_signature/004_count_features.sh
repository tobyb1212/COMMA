#!/bin/bash

while getopts d:o: option
do
	case "${option}" in
		d) DataDir=${OPTARG};;
		o) OutDir=${OPTARG};;
		*) echo "usage: $0 [-d] [-o]" >&2
			exit 1;;
	esac
done

if [ -d "$DataDir" ]
then
	echo "Data directory $DataDir found"
else
	echo "Data directory not found"
	exit 0
fi

if [ -d "$OutDir" ]
then
        echo "OutDir $OutDir found"
else
        echo "Output directory not provided"
        exit 0
fi

cd $DataDir

while IFS=',' read -r srr sample; do
    infile="${srr}Aligned_filtered.sortedByCoord.out.bam"
    outfile="mapped_bam/${sample}.bam"
    if [[ -f "$infile" ]]; then
        echo "Copying $infile → $outfile"
        cp "$infile" "$outfile"
    else
        echo "WARNING: $infile not found!"
    fi
done < /COMMA/results/bulk_rna/pdgf/mapping/mapping.csv

# using featureCounts in Rsubread
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2 

Rscript /COMMA/002_VSMC_phenotypes/004_PDGF_signature/R_scripts/004_Feature_Counts.r -d "$DataDir" -o "$OutDir" 
