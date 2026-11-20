import numpy as np
import pandas as pd
import os
from scipy.stats import norm

os.chdir("/COMMA/results/GWAS/CMVD/finemapping")

SUMSTATS = "/COMMA/results/GWAS/CMVD/sumstats/NTCV_Supplementary_Table4.csv"
BFILE = "/COMMA/reference_genome/ldsc/plink_map/merged/merged_ALL_chr"

THRESHOLD = 1e-5
WINDOW = 500_000

gwas = pd.read_csv(SUMSTATS, sep=r"\s+")
bim = pd.read_csv(BFILE + ".bim", sep=r"\s+", header=None,
                  names=["CHR", "SNP", "CM", "BP", "REF_A1", "REF_A2"])

# Standardise chromosomes and positions
for df in [gwas, bim]:
    df["CHR"] = (df["CHR"].astype(str).str.replace("chr", "", case=False, regex=False)
                 .str.upper().replace({"X": "23"}).astype(int))
    df["BP"] = df["BP"].astype(int)

gwas["SNP"] = gwas["CHR"].astype(str) + ":" + gwas["BP"].astype(str)
bim = bim.rename(columns={"SNP": "BIM_SNP"})

# Match GWAS to BIM
gwas = gwas.merge(bim[["BIM_SNP", "CHR", "BP", "REF_A1", "REF_A2"]],
                  on=["CHR", "BP"], how="inner")

# Harmonise alleles
for col in ["A1", "A2", "REF_A1", "REF_A2"]:
    gwas[col] = gwas[col].str.upper()

gwas = gwas[~(gwas["A1"] + gwas["A2"]).isin(["AT", "TA", "CG", "GC"])].copy()

same = (gwas["A1"] == gwas["REF_A1"]) & (gwas["A2"] == gwas["REF_A2"])
flip = (gwas["A1"] == gwas["REF_A2"]) & (gwas["A2"] == gwas["REF_A1"])
gwas = gwas[same | flip].copy()

# Calculate and align effects
gwas["BETA"] = np.log(gwas["OR"])
gwas["Z"] = norm.isf(np.clip(gwas["P"], np.finfo(float).tiny, 1) / 2) * np.sign(gwas["BETA"])
flip = gwas["A1"] == gwas["REF_A2"]
gwas["BETA_ALIGNED"] = np.where(flip, -gwas["BETA"], gwas["BETA"])
gwas["Z_ALIGNED"] = np.where(flip, -gwas["Z"], gwas["Z"])
gwas["SE"] = np.where(gwas["Z"] != 0, abs(gwas["BETA"] / gwas["Z"]), np.nan)

gwas = (gwas.replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["P", "BETA_ALIGNED", "Z_ALIGNED", "SE"])
        .sort_values("P").drop_duplicates("BIM_SNP"))

gwas.to_csv("genomewide_harmonised_sumstats.tsv.gz",
            sep="\t", index=False, compression="gzip")

# Extract suggestive SNPs after harmonisation
hits = gwas[gwas["P"] < THRESHOLD].sort_values(["CHR", "BP"]).copy()
hits.to_csv("suggestive_snps.tsv", sep="\t", index=False)

# Merge overlapping +/- 500 kb windows
regions = []
for chrom, group in hits.groupby("CHR"):
    start = end = None
    members = []
    for _, snp in group.iterrows():
        s, e = max(1, int(snp["BP"]) - WINDOW), int(snp["BP"]) + WINDOW
        if start is not None and s > end:
            lead = min(members, key=lambda x: x["P"])
            regions.append([chrom, start, end, lead["SNP"],
                            lead["BP"], lead["P"], len(members)])
            members = []
        start = s if not members else start
        end = max(end, e) if members else e
        members.append(snp)
    if members:
        lead = min(members, key=lambda x: x["P"])
        regions.append([chrom, start, end, lead["SNP"],
                        lead["BP"], lead["P"], len(members)])

regions = pd.DataFrame(regions, columns=[
    "CHR", "START", "END", "LEAD_SNP", "LEAD_BP", "LEAD_P", "N_HITS"
])
regions.insert(0, "LOCUS", [f"locus_{i+1:03d}" for i in range(len(regions))])
regions.to_csv("suggestive_loci.tsv", sep="\t", index=False)

print("Suggestive SNPs:", len(hits))
print("Candidate loci:", len(regions))
print("Harmonised SNPs:", len(gwas))