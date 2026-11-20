#!/bin/bash

# Download 1000 Genomes from Phase 3 in hg38
cd /COMMA/reference_genome/ldsc/vcf
BASE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/"
wget -r -np -nH -nd -A "*.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz, *.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi" "$BASE_URL"

export PATH="$PATH:software/plink/plink"
export PATH=software/bcftools-1.19/bin:$PATH

# Download dbSNP157 Files with Accurate Allele Frequencies & rsids
cd /COMMA/data/GWAS/CAD_ANGINA/snps/all/
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi
module load HTSlib/1.18-GCC-12.3.0
( zcat GCF_000001405.40.gz | grep '^#' ;
  zcat GCF_000001405.40.gz | grep -v '^#' | grep -E '(^|;)COMMON($|;)' ) \
  | bgzip -c > dbsnp_common_hg38.vcf.gz
bcftools view -i 'INFO/COMMON' GCF_000001405.40.gz -Oz -o dbsnp_common_hg38.vcf.gz
bcftools annotate --rename-chrs rename_chr.txt -Oz -o dbsnp_common_hg38_fixed.vcf.gz dbsnp_common_hg38.vcf.gz
bcftools index -t dbsnp_common_hg38_fixed.vcf.gz

# Subset 1000 Genomes File to European Ancestry
cd /COMMA/reference_genome/ldsc/vcf
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
awk '($7=="CEU"||$7=="TSI"||$7=="GBR"||$7=="FIN"||$7=="IBS"){print $2}' 20130606_g1k.ped > eur_samples.txt
for chr in {1..22} X; do # Improve This? If possible
  bcftools view -S eur_samples.txt -Oz -o EUR.chr${chr}.vcf.gz \
    ALL.chr${chr}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz --force-samples
  bcftools index -t EUR.chr${chr}.vcf.gz
done

# Download UKBB LDSC Files
wget https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz
gunzip UKBB.ALL.ldscore.tar.gz
tar -xvf UKBB.ALL.ldscore.tar

# Construct Plink Formats for 1000G
for file in EUR*.vcf.gz; do
  base="${file%.vcf.gz}"
  software/plink/plink --vcf "$file" --make-bed --out "../plink/$base"
done

cd /COMMA/reference_genome/ldsc/maps

wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip plink.GRCh38.map.zip 

# Convert map to correct SHAPEIT FORMAT
for file in *.map; do
    {
        echo "pposition rrate gposition"
        awk '{print $4, 0, $3}' "$file"
    } > "${file%.map}.fixed.map"
done

cd /COMMA/reference_genome/ldsc/

# Iterate over every chromosome file and add the cM to it
for chr in {1..22} X; do
    base="EUR.chr${chr}"
    mapfile="maps/plink.chr${chr}.GRCh38.fixed.map"
    out="${base}_map"
    if [[ -f "plink/${base}.bim" && -f "$mapfile" ]]; then
        echo "Processing chromosome ${chr}..."
        software/plink/plink --bfile "plink/$base" \
              --cm-map "$mapfile" $chr --make-bed \
              --out "plink_map/${out}" --set-missing-var-ids @:#
    else
        echo "Skipping chromosome ${chr}: missing .bim or map file"
    fi
done

# Then iterate over each chromosome and perform ldsc on it 
module load Miniforge3
eval "$(conda shell.bash hook)"
conda activate ldsc

for chr in {1..22} X; do
    base="plink_map/EUR.chr${chr}_map"
    out="L2/${chr}_L2"
    if [[ -f "${base}.bim" ]]; then
        echo "Processing chromosome ${chr}..."
    python conda/skylake/envs/ldsc/ldsc/ldsc.py \
    --bfile "$base" --l2 --ld-wind-kb 1000 --out "$out"
    else
        echo "Skipping chromosome ${chr}: missing .bim"
    fi
done

# Merge L2 Output Files (adding the M files together)
zcat L2/1_L2.l2.ldscore.gz | head -n 1 > merged_L2/merged.L2.ldscore
for chr in {1..22}; do
    zcat L2/${chr}_L2.l2.ldscore.gz | tail -n +2 >> merged_L2/merged.L2.ldscore
done
zcat L2/X_L2.l2.ldscore.gz | tail -n +2 >> merged_L2/merged.L2.ldscore
gzip merged_L2/merged.L2.ldscore

# Sum M Values
total_M=0
for chr in {1..22} X; do
    val=$(cat L2/${chr}_L2.l2.M)
    total_M=$(echo "$total_M + $val" | bc)
done
echo $total_M > merged_L2/merged.l2.M

# Sum M_5_50 values
total_M50=0
for chr in {1..22} X; do
    val=$(cat L2/${chr}_L2.l2.M_5_50)
    total_M50=$(echo "$total_M50 + $val" | bc)
done
echo $total_M50 > merged_L2/merged.l2.M_5_50

# Make Simple Annotation File with Merged Bim File and Bed Output
cd ./plink_map
ls *bim | sed 's/.bim//' | tail -n +2 > merge_list.txt
software/plink/plink --bfile EUR.chr1_map --merge-list merge_list.txt --make-bed --out /COMMA/reference_genome/ldsc/plink_map/merged/merged_ALL_chr
cd ..

# Finalise Plink Output
software/plink/plink --bfile /COMMA/reference_genome/ldsc/plink_map/merged/merged_ALL_chr \
  --make-bed --out /COMMA/reference_genome/ldsc/EURhg38/EURhg38
software/plink/plink --bfile /COMMA/reference_genome/ldsc/plink_map/merged/merged_ALL_chr \
  --freq --out /COMMA/reference_genome/ldsc/EURhg38/EURhg38

# Finalise posid Format
cd /COMMA/reference_genome/ldsc/EURhg38
awk 'BEGIN{OFS="\t"}{print $1, $1":"$4, $3, $4, $5, $6}' EURhg38.bim > EUR_posid_hg38.bim
cp EURhg38.bed EUR_posid_hg38.bed
cp EURhg38.fam EUR_posid_hg38.fam
software/plink/plink --bfile EUR_posid_hg38 --freq --out EUR_posid_hg38

# For mtCojo Input - Remove Rare SNPs & Poorly Genotyped SNPs
cd /COMMA/reference_genome/ldsc/vcf
for file in EUR*.vcf.gz; do
  base="${file%.vcf.gz}"
  bcftools annotate -x ID -Oz -o "../mtcojo/id_vcf/${base}.id.vcf.gz" \
    --set-id '%CHROM:%POS:%REF:%ALT' "$file"
  bcftools index "../mtcojo/id_vcf/${base}.id.vcf.gz"
  software/plink/plink --vcf "../mtcojo/id_vcf/${base}.id.vcf.gz" --make-bed --out "../mtcojo/plink/$base" --maf 0.01 --geno 0.05
done
cd /COMMA/reference_genome/ldsc
for chr in {1..22} X; do
    base="EUR.chr${chr}"
    mapfile="maps/plink.chr${chr}.GRCh38.fixed.map"
    out="${base}_map"
    if [[ -f "mtcojo/plink/${base}.bim" && -f "$mapfile" ]]; then
        echo "Processing chromosome ${chr}..."
        software/plink/plink --bfile "mtcojo/plink/${base}" \
              --cm-map "$mapfile" $chr --make-bed \
              --out "mtcojo/plink_map/${out}" --set-missing-var-ids @:#
    else
        echo "Skipping chromosome ${chr}: missing .bim or map file"
    fi
done
for chr in {1..22} X; do
    base="mtcojo/plink_map/EUR.chr${chr}_map"
    out="mtcojo/L2/${chr}"
    if [[ -f "${base}.bim" ]]; then
        echo "Processing chromosome ${chr}..."
    python conda/skylake/envs/ldsc/ldsc/ldsc.py \
    --bfile "$base" --l2 --ld-wind-kb 1000 --out "$out"
    else
        echo "Skipping chromosome ${chr}: missing .bim"
    fi
done
