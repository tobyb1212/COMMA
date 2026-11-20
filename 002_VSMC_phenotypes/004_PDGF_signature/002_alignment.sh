#!/bin/bash

if [ -z "$SLURM_ARRAY_TASK_ID" ]
then    
                echo "Not an array job - need to submit the job as an array in the form of sbatch --array=1-[N_Samples] FastQC_Script_SLURM_1.sh"
                exit 0 
else 
                i=$((SLURM_ARRAY_TASK_ID-1)) 
fi 

while getopts d:o:r: option
do
        case "${option}" in
                d) DataDir=${OPTARG};;
                o) OutDir=${OPTARG};;
                r) RefDir=${OPTARG};;
                *) echo "usage: $0 [-d DataDir] [-o OutputDIR] [-r STAR_Reference]" >&2
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

if [ -d "$RefDir" ] #Reference Directory is the one containing the STAR_index genome.
then
	echo "Reference directory $RefDir found"
else
	echo "Reference directory not found"
	exit 0
fi

cd "$DataDir" || exit

fnames=()
for fname in *_1.fastq.gz
do
	fname=${fname%_1.fastq.gz}
	fnames[${#fnames[@]}]=$fname #Storing the list of samples as a variable
done

# start alignment

seqfile1=$DataDir/${fnames[$i]}_1.fastq.gz #Setting These Two Variables for the Raw Sequences for the Forward and Reverse Mate Pairs
seqfile2=$DataDir/${fnames[$i]}_2.fastq.gz

#check filename
if [ -e "$seqfile1" ] #Checking whether the file exists okay or not.
then
	echo "$seqfile1 found"
else
	echo "$seqfile1 not found"
	exit 0
fi
if [ -e "$seqfile2" ] #Checking whether the file exists okay or not.
then
	echo "$seqfile2 found"
else
	echo "$seqfile2 not found"
	exit 0
fi

#STAR mapping

echo "mapping $seqfile1 using STAR default settings"

#mapping - change paths
outprefix=$OutDir/${fnames[$i]}

module load STAR/2.7.11a-GCC-12.3.0

STAR --runThreadN 2 --genomeDir $RefDir \
	--readFilesIn $seqfile1 $seqfile2 \
	--outFileNamePrefix $outprefix --outSAMunmapped Within \
	--outFilterMultimapNmax 10 --alignSJoverhangMin 5 \
	--alignSJDBoverhangMin 3 --outFilterMismatchNmax 10 \
	--outFilterMismatchNoverLmax 0.04 --alignIntronMin 21 \
	--outFilterIntronMotifs RemoveNoncanonical \
	--outSAMtype BAM SortedByCoordinate --readFilesCommand zcat

module load SAMtools/1.18-GCC-12.3.0

samtools view -h -q 255 -b $OutDir/${fnames[$i]}Aligned.sortedByCoord.out.bam > $OutDir/${fnames[$i]}Aligned_filtered.sortedByCoord.out.bam

rm $OutDir/${fnames[$i]}Aligned.sortedByCoord.out.bam
