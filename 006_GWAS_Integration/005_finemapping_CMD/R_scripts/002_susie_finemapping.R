library(data.table)
library(susieR)

setwd("/COMMA/results/GWAS/CMVD/finemapping")

BFILE <- "/COMMA/reference_genome/ldsc/plink_map/merged/merged_ALL_chr"
PLINK <- "/COMMA/software/plink/plink"

regions <- fread("suggestive_loci.tsv")
gwas <- fread("genomewide_harmonised_sumstats.tsv.gz")
gwas[, N := 1002]

dir.create("susie_results", showWarnings=FALSE)

results <- lapply(seq_len(nrow(regions)), function(i) {
  x <- regions[i]
  id <- x$LOCUS
  if (x$CHR == 23) return(NULL)

  dir.create(paste0("susie_results/", id), showWarnings=FALSE)
  prefix <- paste0("susie_results/", id, "/", id)

  tryCatch({
    sub <- gwas[CHR == x$CHR & BP >= x$START & BP <= x$END &
                is.finite(Z_ALIGNED)]

    if (nrow(sub) < 10 || nrow(sub) > 10000) return(NULL)

    # Extract SNPs using PLINK IDs
    fwrite(sub[, .(BIM_SNP)], paste0(prefix, ".snps"), col.names=FALSE)

    status <- system2(PLINK, c(
      "--bfile", BFILE, "--chr", x$CHR,
      "--from-bp", x$START, "--to-bp", x$END,
      "--extract", paste0(prefix, ".snps"),
      "--maf", "0.001", "--geno", "0.05",
      "--make-bed", "--out", prefix
    ))
    if (status != 0) stop("PLINK extraction failed")

    # Match GWAS SNP order to PLINK
    bim <- fread(paste0(prefix, ".bim"), header=FALSE)
    sub <- sub[match(bim$V2, sub$BIM_SNP)]

    stopifnot(
      !anyNA(sub$BIM_SNP),
      all(sub$BIM_SNP == bim$V2)
    )

    # Align Z-scores to extracted PLINK A1
    same <- sub$REF_A1 == bim$V5 & sub$REF_A2 == bim$V6
    flip <- sub$REF_A1 == bim$V6 & sub$REF_A2 == bim$V5
    stopifnot(all(same | flip))

    z <- ifelse(flip, -sub$Z_ALIGNED, sub$Z_ALIGNED)

    # Calculate LD matrix
    status <- system2(PLINK, c(
      "--bfile", prefix, "--r", "square",
      "--out", prefix
    ))
    if (status != 0) stop("PLINK LD calculation failed")

    R <- as.matrix(fread(paste0(prefix, ".ld"), header=FALSE))
    stopifnot(nrow(R) == length(z), ncol(R) == length(z))

    # Run SuSiE
    fit <- susie_rss(z=z, R=R, n=1002, L=10)

    saveRDS(fit, paste0(prefix, ".susie.rds"))

    fwrite(data.table(
      SNP=sub$BIM_SNP,
      PIP=fit$pip
    ), paste0(prefix, ".pip.tsv"), sep="\t")

    fwrite(data.table(
      LOCUS=id,
      N_SNPS=length(z),
      N_CS=length(fit$sets$cs),
      CONVERGED=fit$converged
    ), paste0(prefix, ".summary.tsv"), sep="\t")

    message(id, ": completed")
    return(fit)

  }, error=function(e) {
    message(id, ": ", conditionMessage(e))
    return(NULL)
  })
})

cat("Completed loci:", sum(!sapply(results, is.null)), "/", nrow(regions), "\n")

# Extract Lead Loci
files <- list.files("susie_results", pattern="\\.pip.tsv$",
                    recursive=TRUE, full.names=TRUE)

leads <- rbindlist(lapply(files, function(f) {
  x <- fread(f)
  x <- x[which.max(PIP)]
  x[, LOCUS := basename(dirname(f))]
  x
}), fill=TRUE)

setcolorder(leads, c("LOCUS", "SNP", "PIP"))
fwrite(leads, "susie_lead_variants.tsv", sep="\t")

print(leads)