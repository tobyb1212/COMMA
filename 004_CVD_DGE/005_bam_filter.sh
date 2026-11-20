#!/bin/bash

module load Miniforge3
eval "$(conda shell.bash hook)"
conda activate scanpy 

python /COMMA/004_CVD_DGE/Python_scripts/005_barcode_filter.py

if [ -z "$SLURM_ARRAY_TASK_ID" ]
then 
                echo "Not an array job! This script needs to be submitted as an array job with sbatch --array=1-[N_samples] "
                exit 0
else 
                i=$((SLURM_ARRAY_TASK_ID-1))
fi

module load SAMtools
module load deepTools

cd //COMMA/data/processed/rna_koenig/bam_filtering

mv reruns/*.bam .

fnames=()
for fname in *.txt
do
fnames[${#fnames[@]}]=${fname%.txt}
done

samtools view -b -D CB:"${fnames[$i]}.txt" "${fnames[$i]}Aligned.sortedByCoord.out.bam" | samtools sort -o "${fnames[$i]}.filtered.bam"
samtools index "${fnames[$i]}.filtered.bam"
bamCoverage -b "${fnames[$i]}.filtered.bam" -o "${fnames[$i]}.bw" --normalizeUsing CPM -p 6
