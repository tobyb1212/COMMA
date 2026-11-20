#!/bin/bash

FOLDER="/COMMA/results/scrna/downstream_analysis/epicardial_coronary/data/"
cd $FOLDER

mkdir $FOLDER/concatenated

prefixes=$(ls ./*.fastq.gz 2>/dev/null | sed -E 's/_L[0-9]+_R[12]_[0-9]+\.fastq\.gz$//' | sed -E 's/_R[12]_[0-9]+\.fastq\.gz$//' | sed -E 's/\.fastq\.gz$//' | sort | uniq)

for prefix in $prefixes; do
    r1_files=($(ls "${prefix}"*_R1_*.fastq.gz 2>/dev/null | sort -V))
    r2_files=($(ls "${prefix}"*_R2_*.fastq.gz 2>/dev/null | sort -V))
    i1_files=($(ls "${prefix}"*_I1_*.fastq.gz 2>/dev/null | sort -V))
    output_r1="${FOLDER}/concatenated/$(basename ${prefix})_L001_R1_001.fastq.gz"
    output_r2="${FOLDER}/concatenated/$(basename ${prefix})_L001_R2_001.fastq.gz"
    output_i1="${FOLDER}/concatenated/$(basename ${prefix})_L001_I1_001.fastq.gz"
    if [ ${#r1_files[@]} -gt 1 ]; then
        cat "${r1_files[@]}" > "$output_r1"
        echo "Concatenated R1 files for ${prefix} into $output_r1"
    elif [ ${#r1_files[@]} -eq 1 ]; then
        cp "${r1_files[0]}" "$output_r1"
        echo "Moved R1 file for ${prefix} to $output_r1"
    fi
    if [ ${#r2_files[@]} -gt 1 ]; then
        cat "${r2_files[@]}" > "$output_r2"
        echo "Concatenated R2 files for ${prefix} into $output_r2"
    elif [ ${#r2_files[@]} -eq 1 ]; then
        cp "${r2_files[0]}" "$output_r2"
        echo "Moved R2 file for ${prefix} to $output_r2"
    fi
    if [ ${#i1_files[@]} -gt 1 ]; then
        cat "${i1_files[@]}" > "$output_i1"
        echo "Concatenated I1 files for ${prefix} into $output_i1"
    elif [ ${#i1_files[@]} -eq 1 ]; then
        cp "${i1_files[0]}" "$output_i1"
        echo "Moved I1 file for ${prefix} to $output_i1"
    fi
done

