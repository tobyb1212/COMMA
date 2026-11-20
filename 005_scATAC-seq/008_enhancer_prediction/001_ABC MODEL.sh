#!/bin/bash

cd /COMMA/results/scatac/multiomic_integration/scE2G/vascular/fragments

# Create and Format Pseudobulk Fragment File
LC_ALL=C zcat 'Vascular Smooth Muscle Cell.tsv.gz' | sed '/^#/d' | awk -v OFS='\t' '{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+"; print $1,mid+1,$3,"N",1000,"-"}' | sort -k 1,1V -k 2,2n -k3,3n --parallel 5 | bgzip -c > VSMC.tagAlign.gz  
tabix -p bed VSMC.tagAlign.gz
cp VSMC.tagAlign.gz /COMMA/results/scatac/multiomic_integration/ABC/vsmc

# Compress ChIP-seq Alignments Normalised by Library Size
cd /COMMA/results/chIP/bam/deduplicated
MIN=$(for f in H3K27ac_*.rm_dup.bam; do samtools view -c -f 2 -F 1804 "$f"; done | sort -n | head -1)
for f in H3K27ac_*.rm_dup.bam; do N=$(samtools view -c -f 2 -F 1804 "$f"); F=$(awk -v m=$MIN -v n=$N 'BEGIN{printf "%.6f",m/n}'); samtools view -bs 42${F#0} "$f" > "${f%.bam}.norm.bam"; done
samtools merge -@ 6 H3K27ac_merged.bam H3K27ac_*.norm.bam
samtools sort -@ 6 -o H3K27ac_merged.sorted.bam H3K27ac_merged.bam
samtools index H3K27ac_merged.sorted.bam

# Build Hi-C .hic File with Juicer
cd /COMMA/results/CCC/CA-VSMCs/results
pairtools parse --chroms-path GRCh38.chrom.sizes HiCUP_filtered.bam \
  | pairtools sort \
  | pairtools dedup \
  | pairtools select '(pair_type=="UU")' \
  > contact_map.pairs
grep -v '^#' contact_map.pairs | \
awk 'BEGIN{OFS="\t"} {print 0,$2,$3,0,1,$4,$5,1}' > contact_map.juicer.txt
java -Xmx32g -jar software/juicer/juicer_tools_1.22.01.jar pre \
  -j 6 contact_map.juicer.txt Juicer/contact_map.hic GRCh38.chrom.sizes

# Supplemental Table 11
cd software/ABC-Enhancer-Gene-Prediction
snakemake -j1