#!/bin/bash

if [ -z "$SLURM_ARRAY_TASK_ID" ]
then 
                echo "Not an array job! This script needs to be submitted as an array job with sbatch --array=1-[N_samples] "
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
                *) echo "usage: $0 [-d DataDir] [-o OutDIR] [-r RefDir]" >&2
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

if [ -d "$RefDir" ] #Reference Directory is the one containing the Bowtie_index files.
then
	echo "Reference directory $RefDir found"
else
	echo "Reference directory not found"
	exit 0
fi

cd "$DataDir" || exit

fnames=()
for fname in *_1.fq.gz
do
	fname=${fname%_1.fq.gz}
	fnames[${#fnames[@]}]=$fname #Storing the list of samples as a variable
done

# start alignment

seqfile1=$DataDir/${fnames[$i]}_1.fq.gz #Setting These Two Variables for the Raw Sequences for the Forward and Reverse Mate Pairs
seqfile2=$DataDir/${fnames[$i]}_2.fq.gz

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

echo "mapping $seqfile1 and $seqfile2 using Bowtie2 default settings"

#mapping - change paths
outprefix=$OutDir/${fnames[$i]}

module load Bowtie2
module load SAMtools

bowtie2 --very-sensitive -x "$RefDir"/hg38 -X 2000 -1 "$seqfile1" -2 "$seqfile2" -p 4 \
  |  samtools view -u - \
  |  samtools sort -o "${outprefix}.sortbyC.bam"  -

samtools index -@ 4 "${outprefix}.sortbyC.bam"
