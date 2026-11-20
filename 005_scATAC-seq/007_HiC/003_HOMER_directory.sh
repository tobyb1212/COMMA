#!/bin/bash

i=$((SLURM_ARRAY_TASK_ID-1))

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

cd "$DataDir" || exit

fnames=()
for fname in *.hicup.bam
do
fnames[${#fnames[@]}]=${fname%.hicup.bam}
done

export PATH="$PATH:software/HiCUP-0.9.2/Conversion"
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
module load SAMtools/1.18-GCC-12.3.0

hicup2homer $DataDir/${fnames[$i]}.hicup.bam

export PATH="$PATH:software/HOMER/bin"
module load SAMtools/1.18-GCC-12.3.0

cd "$OutDir" || exit
echo "processing ${fnames[$i]}......"

makeTagDirectory ${fnames[$i]} -format HiCsummary $DataDir/${fnames[$i]}.hicup.bam.homer