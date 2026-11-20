#!/bin/bash

#  sbatch --array=1-12  --export=DataDir=/COMMA/data/meta_analysis/fastq/atac_selewa/concatenated,OutDir=/COMMA/data/meta_analysis/processed/atac_selewa   /COMMA/001_alignment/004_cellranger_atac_secondary.sh

export PATH="$PATH:software/cellranger-atac-2.1.0/bin"

i=$((SLURM_ARRAY_TASK_ID-1))

while getopts ":d:o:" option; do 
    case "$option" in 
        d) DataDir=${OPTARG} ;;
        o) OutDir=${OPTARG} ;;
        *) echo "Usage: $0 [-d input_directory] [-o output_directory]" >&2 
           exit 1 ;;
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
        echo "Output directory $OutDir found"
else
        echo "Output directory not found"
        exit 0
fi

cd "$DataDir" || exit

fnames=()
for fname in *_S1_L001_R1_001.fastq.gz
do
fnames[${#fnames[@]}]=${fname%_S1_L001_R1_001.fastq.gz}
done

id_val="${fnames[$i]:0:64}"

# Using CellRanger-ATAC for all 10X ATAC-seq only Samples

cd "$OutDir" || exit
echo "processing ${fnames[$i]}......"
if ! cellranger-atac count --id "$id_val" \
    --reference /COMMA/reference_genome/cellranger_meta/cellranger-arc_gencode46_jaspar2024/GRCh38/ \
    --sample ${fnames[$i]} --fastqs $DataDir --localcores=8; then
    # Same Reference as for CellRanger-arc 
    echo "${fnames[$i]}" >> "$DataDir/failed_cellranger.txt" # Note if CellRanger-ATAC fails
    exit 1
fi