library(tidyverse)
library(optparse)
library(data.table)

source("/COMMA/006_GWAS_Integration/005_finemapping_CMD/R_scripts/GWAS_FUNCTIONS.R")

plink = "/COMMA/software/plink/plink"
vcf = "/COMMA/data/GWAS/CAD_ANGINA/snps/all/dbsnp_common_hg38_fixed.vcf.gz"
genotypeData = "/COMMA/reference_genome/ldsc/EURhg38/EUR_posid_hg38"

cmvd <- read.csv("/COMMA/results/GWAS/CMVD/finemapping/susie_lead_variants.tsv", sep = "\t")
cmvd_susie <- cmvd$SNP
gwas.hg38.finalist <- getld_GWAS_SNP(cmvd_susie, plink = plink, genotypeData = genotypeData, r2 = 0.8, return_clump=F)
gwas.hg38.finalist <- pos2rsid(gwas.hg38.finalist, bcftools="bcftools", vcf=vcf)
gwas.hg38.gr <- snp2gr(gwas.hg38.finalist)

peak_files <- list.files("/COMMA/results/meta_analysis/scatac/gwas_overlaps/vascular",pattern = "\\.bed$",full.names = TRUE)

overlaps <- do.call(rbind,
    Filter(Negate(is.null), lapply(peak_files, function(pf){
      print(pf)
      peaks <- read.table(pf, sep = "\t")
      genomicstring_peaks <- paste(peaks$V1, peaks$V2, peaks$V3, sep = "-")
      res <- try(GWASOverlaps(SNPList = gwas.hg38.gr, peakSet = genomicstring_peaks,
        setName = tools::file_path_sans_ext(basename(pf))), silent = TRUE)
      if (inherits(res, "try-error") || is.null(res) || NROW(res) == 0) return(NULL)
      res}))) 
write.table(overlaps, "cmvd_gwas_overlaps.tsv", sep = "\t", quote = FALSE, row.names = FALSE)