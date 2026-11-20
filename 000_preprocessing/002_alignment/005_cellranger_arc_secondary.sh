#!/bin/bash

export PATH="$PATH:software/cellranger-arc-2.0.2/bin"

i=$((SLURM_ARRAY_TASK_ID-1))

while getopts "d:o:" option; do 
    case "$option" in 
        d) LibraryDir=${OPTARG} ;;
        o) OutDir=${OPTARG} ;;
        *) echo "Usage: $0 [-d LibraryDir] [-o OutDir]" >&2
            exit 1
            ;;
    esac
done

if [ -d "$LibraryDir" ]
then
        echo "Data Directory $LibraryDir found" 
else
        echo "Data Directory not found"
        exit 0 
fi 

if [ -d "$OutDir" ] 
then 
        echo "OutDirectory $OutDir found" 
else
        echo "OutDirectory being created"
        mkdir"$OutDir" 
fi 

cd $LibraryDir

lnames=()
for lname in *.csv
do
lnames[${#lnames[@]}]=${lname%.csv}
done

# Libraries are .csv files in the following format
# fastqs,sample,library_type
# File Path,Sample,Gene Expression
# File Path,Sample,Chromatin Accessibility
# This is a necessary input to CellRanger-arc

echo ${lnames[$i]}

i=$((SLURM_ARRAY_TASK_ID-1)) 

cd "$OutDir" || exit
echo "processing ${lnames[$i]}"
if ! cellranger-arc count --id ${lnames[$i]} --reference '/COMMA/reference_genome/cellranger_meta/cellranger-arc_gencode46_jaspar2024/GRCh38'  \
    --libraries $LibraryDir/${lnames[$i]}.csv --no-bam --localcores=8; then
    # Same Reference as for CellRanger-atac 
    echo "${fnames[$i]}" >> "$DataDir/failed_cellranger.txt" # Note if CellRanger-arc fails
    exit 1
fi