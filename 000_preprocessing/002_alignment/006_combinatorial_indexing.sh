#!/bin/bash

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

module load Miniforge3
eval "$(conda shell.bash hook)"
conda activate interstellar

cd "$DataDir" || exit

fnames=()
for fname in *_S1_L001_R1_001.fastq.gz
do
fnames[${#fnames[@]}]=${fname%_S1_L001_R1_001.fastq.gz}
done

python /COMMA/000_preprocessing/002_alignment/Python_scripts/006_combinatorial_indexing.py \
  -i ${fnames[$i]}_S1_L001_R2_001.fastq.gz -d $DataDir \
  -b software/cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-cratac-v1.txt \
  -o ${fnames[$i]} # Create I2 fastq file for ${fnames[$i]} using new barcodes that are compatible with CellRanger atac

echo "Converted ${fnames[$i]}....."
conda deactivate

export PATH="$PATH:software/cellranger-atac-2.1.0/bin"

cd "$OutDir" || exit

if ! cellranger-atac count --id ${fnames[$i]} \
   --reference /COMMA/reference_genome/cellranger_meta/cellranger-arc_gencode46_jaspar2024/GRCh38/ \
   --sample ${fnames[$i]} --fastqs $DataDir --localcores=8; then
      # Same Reference as for CellRanger-arc 
      echo "${fnames[$i]}" >> "$DataDir/failed_cellranger.txt" # Note if CellRanger-ATAC fails
      exit 1
fi
