#!/bin/bash

if [ -z "$SLURM_ARRAY_TASK_ID" ]
then    
                echo "Not an array job - need to submit the job as an array in the form of sbatch --array=1-[N_Samples] FastQC_Script_SLURM_1.sh"
                exit 0 
else 
                i=$((SLURM_ARRAY_TASK_ID-1)) 
fi 

while getopts d:o: option # Defining the Input & Output Directories Using GetOpts
do 
        case "{$option}" in 
                d) DataDir=${OPTARG};;
                o) OutDir=${OPTARG};;
                *) echo "usage: $0 [-d] [-o]" >&2 
                        exit 1;; 
        esac
done

if [ -d "$DataDir" ] 
then
        echo "Data Directory $DataDir found" 
else
        echo "Data Directory not found"
        exit 0 
fi 

if [ -d "$OutDir" ] 
then 
        echo "OutDirectory $OutDir found" 
else
        echo "OutDirectory being created"
        mkdir "$OutDir"
fi

cd "$DataDir" || exit 
fnames=()
for file in *.fastq.gz
do 
        fnames[${#fnames[@]}]=$file
done
echo "FastQC on ${fnames[$i]}"

module load FastQC/0.11.9-Java-11 

fastqc "${fnames[$i]}" -o "$OutDir" 