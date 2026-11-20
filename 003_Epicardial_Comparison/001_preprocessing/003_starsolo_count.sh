#!/bin/bash

i=$((SLURM_ARRAY_TASK_ID-1))

module load STAR/2.7.11b-GCC-12.3.0

# Parsing Sequencing Technology from Command Line -v Option

while getopts "d:o:v:" option; do 
    case "$option" in 
        d) DataDir=${OPTARG} ;;
        o) OutDir=${OPTARG} ;;
        v) Tech=${OPTARG} ;;
        *) echo "Usage: $0 [-d input_directory] [-o output_directory] [-v technology]" >&2 
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

if [ -n "$Tech" ]
then
        echo "Techonology $Tech found"
else
        echo "Technology not found"
        exit 0
fi

cd "$DataDir" || exit

fnames=()
for fname in *_L001_R1_001.fastq.gz
do
fnames[${#fnames[@]}]=${fname%_L001_R1_001.fastq.gz}
done

seqfile2=$DataDir/${fnames[$i]}_L00*_R1_001.fastq.gz # StarSolo requires R1 Read as the Second File
if [ $(echo "$seqfile2" | wc -l) -ne 1 ]; then
    echo "Error: Multiple or no files match pattern for seqfile2: ${DataDir}/${fnames[$i]}_L001_R1_001.fastq.gz"
    exit 1
fi
seqfile1=$DataDir/${fnames[$i]}_L00*_R2_001.fastq.gz
if [ $(echo "$seqfile1" | wc -l) -ne 1 ]; then
    echo "Error: Multiple or no files match pattern for seqfile1: ${DataDir}/${fnames[$i]}_L001_R2_001.fastq.gz"
    exit 1
fi
UMILength=$(if [[ "$Tech" == "3 V2" || "$Tech" == "5 V1" ]]; then echo 10; elif [[ "$Tech" == "Multiome" || "$Tech" == "3 V3" || "$Tech" == "Visium V1" ]]; then echo 12; fi) # UMI Length is Dependent on Sequencing Technology
CBWhitelist=$(if [[ "$Tech" == "3 V2" || "$Tech" == "5 V1" ]]; then 
        echo "software/cellranger-8.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt"
    elif [[ "$Tech" == "Multiome" ]]; then echo "software/cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt"
    elif [[ "$Tech" == "3 V3" ]]; then echo "software/cellranger-8.0.1/lib/python/cellranger/barcodes/3M-february-2018_star.txt"
    elif [[ "$Tech" == "Visium V1" ]]; then echo "software/spaceranger-3.1.2/lib/python/cellranger/barcodes/visium-v1.txt"; 
    else echo "Error: Unsupported Tech value '$Tech'"; exit 1; fi) # Barcode List varies for Technology
BarcodeLength=16 # Barcode Length is constant for the technologies used in this meta-analysis
outprefix=$OutDir/${fnames[$i]}

cd "$OutDir" || exit
echo "processing ${fnames[$i]}......"

if [[ $Tech == "5 V1" ]]
then
        echo "STARsolo running in paired-end mode"
  if ! STAR --runThreadN 6 --genomeDir /COMMA/reference_genome/STAR/STAR_meta --readFilesIn $seqfile2 $seqfile1 --clip5pNbases 39 0 \
    --soloType CB_UMI_Simple --soloBarcodeMate 1 --soloCBwhitelist $CBWhitelist --soloUMIlen $UMILength --outSAMtype None \
    --readFilesCommand zcat --outFilterScoreMin 30 --outFileNamePrefix $outprefix  --soloFeatures Gene GeneFull GeneFull_Ex50pAS --runDirPerm All_RWX \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloMultiMappers EM \
    --soloCBlen $BarcodeLength --soloUMIstart $((BarcodeLength+1)) --soloCellFilter EmptyDrops_CR --soloStrand Forward; then
     echo "${fnames[$i]}" >> "$DataDir/failed_starsolo.txt"
    exit 1
    # StarSolo settings adapted from https://github.com/cellgeni/STARsolo
    # Velocyto for 5 PE Technology is Performed Separately
    fi
else
        echo "STARsolo running in single-end mode for Velocyto output"
        if ! STAR --runThreadN 6 --genomeDir /COMMA/reference_genome/STAR/STAR_meta --readFilesIn $seqfile1 $seqfile2 \
    --soloType CB_UMI_Simple  --soloCBwhitelist $CBWhitelist --outSAMtype None --soloUMIlen $UMILength --runDirPerm All_RWX \
    --readFilesCommand zcat --clipAdapterType CellRanger4 --outFilterScoreMin 30 --outFileNamePrefix $outprefix  --soloFeatures Gene GeneFull Velocyto GeneFull_Ex50pAS \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloMultiMappers EM \
    --soloCBlen $BarcodeLength --soloUMIstart $((BarcodeLength+1)) --soloCellFilter EmptyDrops_CR --soloStrand Forward --soloBarcodeReadLength 0; then
    echo "${fnames[$i]}" >> "$DataDir/failed_starsolo.txt"
    exit 1
    # StarSolo Settings adapted from https://github.com/cellgeni/STARsolo 
    fi
fi
