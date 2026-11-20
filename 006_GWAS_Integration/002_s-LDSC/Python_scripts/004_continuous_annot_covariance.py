import gzip
import numpy as np
import pandas as pd
from pathlib import Path

base = Path("COMMA")
ref = base / "reference_genome/ldsc/plink_map/merged/merged_ALL_chr"
out = base / "results/GWAS/s-LDSC/partitioned_heritability_all"

names = ["Arterial_Endothelial_Cell", "Vascular_Smooth_Muscle_Cell",
         "Capillary_Endothelial_Cell", "Pericyte", "Venous_Endothelial_Cell", "All_OCR"]

# Match common SNPs to BIM order
bim = pd.read_csv(str(ref)+".bim", sep=r"\s+", header=None, usecols=[1])
frq = pd.read_csv(str(ref)+".frq", sep=r"\s+")

common = set(frq.loc[(frq.MAF >= 0.05) & (frq.MAF <= 0.5), "SNP"])
mask = bim[1].isin(common).to_numpy()

assert mask.sum() == 5698384 # N SNPs in LDSC
print("Common SNPs:", mask.sum())

results = []
X = []

for name in names:
    with gzip.open(out / (name+".annot.gz"), "rt") as f:
        assert next(f).strip() == "ANNOT"
        values = np.fromiter((float(line) for line in f), dtype=float)
    assert len(values) == len(mask)
    x = values[mask]
    expected = float((out / (name+".l2.M_5_50")).read_text().split()[0])
    results.append({
        "Annotation": name,
        "N": len(x),
        "Sum": x.sum(),
        "LDSC_Sum": expected,
        "Difference": x.sum()-expected,
        "Mean": x.mean(),
        "SD": x.std(ddof=0)
    })
    X.append(x)
    print(name, "SD =", round(x.std(ddof=0), 9),
          "Difference =", round(x.sum()-expected, 6))

pd.DataFrame(results).to_csv(out / "annotation_validation.csv", index=False)

cov = np.cov(np.column_stack(X), rowvar=False, ddof=0)
pd.DataFrame(cov, index=names, columns=names).to_csv(
    out / "annotation_covariance.csv"
)

print("Saved annotation_validation.csv and annotation_covariance.csv")