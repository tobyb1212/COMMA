#!/bin/bash

if [ -z "$SLURM_ARRAY_TASK_ID" ]
then 
                echo "Not an array job! This script needs to be submitted as an array job with sbatch --array=1-[N_samples] "
                exit 0
else 
                i=$((SLURM_ARRAY_TASK_ID-1))
fi

while getopts d:o: option
do
        case "${option}" in
                d) DataDir=${OPTARG};;
                o) OutDir=${OPTARG};;
                *) echo "usage: $0 [-d DataDir] [-o OutDIR]" >&2
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
	echo "Output directory $OutDir found"
else
	echo "Output directory not found"
	exit 0
fi

cd "$DataDir" || exit

fnames=()
for fname in *.sortbyC.bam
do
        fname=${fname%.sortbyC.bam}
        fnames[${#fnames[@]}]=$fname
done

cd "$OutputDir" || exit 
module load SAMtools
if [ -e "$DataDir/${fnames[$i]}.sortbyC.bam.bai" ] #Seeing if the Files in the Data Directory Exist
then
        echo "Sorted bam file ${fnames[$i]}.sortbyC.bam found"
        echo "Marking duplicated reads"
        java -jar software/picard/picard.jar MarkDuplicates -I "${fnames[$i]}".sortbyC.bam -O "${fnames[$i]}".sortbyC.mk_dup.bam -M "${fnames[$i]}".mk_dup.txt # This Marks the Duplicates with the BAM flag 1024 indicating Duplication, Also creates a text output file with the metrics
        echo "Removing duplicated reads"
        samtools view -@ 3 -b -F 1024 -f 2 "${fnames[$i]}".sortbyC.mk_dup.bam > "${fnames[$i]}".sortbyC.rm_dup.bam # Using 3 Cores, Outputting in .bam and removing any reads flagged as duplicates.
        rm "${fnames[$i]}".sortbyC.mk_dup.bam # Removing this Intermediate File
        echo "Creating new index"
        samtools index -@ 3 "${fnames[$i]}".sortbyC.rm_dup.bam # Creating New Index & Using Three 
else
        echo "Sorted bam file not found."
        exit 0
fi