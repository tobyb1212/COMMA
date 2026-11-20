#!/bin/bash

# Unharmonised CAD SNPs from European Ancestry (UKBB) for https://www.ebi.ac.uk/gwas/studies/GCST90473536
cd /COMMA/data/GWAS/CAD_ANGINA/snps/cad/ukbb
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90473001-GCST90474000/GCST90473536/GCST90473536.tsv.gz-meta.yaml
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90473001-GCST90474000/GCST90473536/GCST90473536.tsv.gz 
gunzip GCST90473536.tsv.gz

# Unharmonised Angina SNPs from European Ancestry (UKBB) for https://www.ebi.ac.uk/gwas/studies/GCST90473533
cd /COMMA/data/GWAS/CAD_ANGINA/snps/angina/ukbb
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90473001-GCST90474000/GCST90473533/GCST90473533.tsv.gz-meta.yaml
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90473001-GCST90474000/GCST90473533/GCST90473533.tsv.gz
gunzip GCST90473533.tsv.gz

# Construct UKBB Tables - SNPs do not need to be converted as They are in hg38 Format Already
awk 'BEGIN {OFS="\t"} 
#NR==1 {print "SNP","A1","A2","freq","b","se","p","N"; next} 
#{print "chr"$9,$3,$4,$7,$5,$6,$8,$11}' /COMMA/data/GWAS/CAD_ANGINA/snps/angina/ukbb/GCST90473533.tsv \
> /COMMA/results/GWAS/mtCOJO/ukbb_angina_test.raw
awk 'BEGIN {OFS="\t"} 
#NR==1 {print "SNP","A1","A2","freq","b","se","p","N"; next} 
#{print "chr"$9,$3,$4,$7,$5,$6,$8,$11}' /COMMA/data/GWAS/CAD_ANGINA/snps/cad/ukbb/GCST90473536.tsv \
> /COMMA/results/GWAS/mtCOJO/ukbb_cad_test.raw

# Setup List of Summary Statistics
ukbb_cad=/COMMA/results/GWAS/mtCOJO/ukbb_cad_test.raw
ukbb_ang=/COMMA/results/GWAS/mtCOJO/ukbb_angina_test.raw
ukbb_cad_freq=0.050
ukbb_ang_freq=0.076
combos=(
"cadMI_ukbb_vs_ang_ukbb|$ukbb_cad|$ukbb_ang|$ukbb_cad_freq|$ukbb_ang_freq")

cd /COMMA/results/GWAS/mtCOJO
mtcojo_ref=/COMMA/results/GWAS/mtCOJO/mtcojo_ref_data.txt
gsmr_ref=/COMMA/results/GWAS/gsmr/gsmr_ref_data.txt

module load R
module load Miniforge3
eval "$(conda shell.bash hook)"
conda activate gcta-mkl
export MKL_ENABLE_INSTRUCTIONS=AVX2
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

for combo in "${combos[@]}"; do
    IFS='|' read -r name cad_sumstats ang_sumstats cad_freq ang_freq <<< "$combo"
    echo "Processing $name"
    mkdir /COMMA/results/GWAS/gsmr/$name
    mkdir /COMMA/results/GWAS/mtCOJO/$name
    cd /COMMA/results/GWAS/gsmr/$name
    # Create GSMR Files
    printf "cad\t%s\n" "$cad_sumstats" > cad.txt
    printf "ang\t%s\n" "$ang_sumstats" > ang.txt
    # Run GSMR
    software/gcta/gcta-1.95.0-linux-kernel-3-x86_64/gcta64 --mbfile $gsmr_ref \
        --gsmr-file cad.txt ang.txt --gsmr-direction 0 --out ${name}_gsmr_result --effect-plot
    awk 'NR>1 {print $1, $3}' ${name}_gsmr_result.gsmr > ../../mtCOJO/$name/${name}_mtcojo_bxy.txt
    # Call R Code to Create Relevant Plots
    # Supplemental Figure 14A
    Rscript - <<RS
        library(tidyverse)
        source("software/gcta/gcta-1.95.0-linux-kernel-3-x86_64/gsmr_plot.r")
        gsmr_data <- read_gsmr_data("${name}_gsmr_result.eff_plot.gz")
        gsmr_summary(gsmr_data)
        pdf("${name}_gsmr.pdf", width = 8, height = 5)
        plot_gsmr_effect(gsmr_data, "cad", "ang", colors()[75])
        dev.off()
RS
    # Create mtCojo Files
    cd /COMMA/results/GWAS/mtCOJO/$name
    mtcojo_list=mtcojo_summary_data.list
    printf "ang %s\n" "$ang_sumstats"  > "$mtcojo_list"
    printf "cad %s" "$cad_sumstats" >> "$mtcojo_list"
    # Run mtCojo
    software/gcta/gcta-1.95.0-linux-kernel-3-x86_64/gcta64 --mbfile $mtcojo_ref --mtcojo-file mtcojo_summary_data.list --ref-ld-chr /COMMA/reference_genome/ldsc/mtcojo/L2/ \
        --w-ld-chr /COMMA/reference_genome/ldsc/mtcojo/L2/ --out ${name}_mtcojo_result --mtcojo-bxy ${name}_mtcojo_bxy.txt
    cd ..
done