import os
import scanpy as sc

vsmc_adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/vsmc_annotated_final.h5ad")

fastq_names = {
    "TWCM-11-74": "SRR15835879_SRR15835882_SRR15835880_SRR15835881",
    "TWCM-13-96": "SRR15835839_SRR15835840_SRR15835841_SRR15835838",
    "TWCM-13-168": "SRR15835910_SRR15835908_SRR15835909_SRR15835911",
    "H_ZC-LVAD": "SRR15835822_SRR15835823",
    "TWCM-13-47": "SRR15835921",
    "TWCM-13-208": "SRR15835857_SRR15835855_SRR15835856_SRR15835854"
}

for sample, fastq in fastq_names.items():
    mask = (
        (vsmc_adata.obs["Study"] == "Koenig 2022") &
        (vsmc_adata.obs["Sample Name"] == sample) &
        (vsmc_adata.obs["fastq_name"] == fastq)
    )
    barcodes = vsmc_adata.obs.loc[mask].index.astype(str).str.rsplit("_", n=1).str[-1]
    filename = f"/COMMA/data/meta_analysis/processed/rna_koenig/bam_filtering/{fastq}.txt"
    with open(filename, "w") as f:
        f.write("\n".join(barcodes) + "\n" if len(barcodes) else "")
    print(f"{sample}: {len(barcodes)} barcodes -> {filename}")