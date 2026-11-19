# COMMA: COronary Microvascular Multiomic Atlas

This repository contains the scripts necessary to reproduce our integrated analysis of human single-cell and single-nucleus datasets to create COMMA. 

---

## 📑 Table of Contents 

1. [🧬 Overview](#️-overview)
2. [👩‍💻 Code Structure](#-code-structure)
3. [🔁 How to Reproduce](#-how-to-reproduce)
4. [📋 Metadata](#-metadata) 
5. [📚 References](#-references)
6. [📝 License](#-license)
7. [💻 Contact](#-contact)

---

## 🧬 Overview

In this repository, we include code to generate the sc/snRNA-seq and scATAC-seq components of COMMA. 

Additionally, code is available for our fine-grained subclustering of microvascular VSMC phenotypes and the comparison of microvascular and epicardial VSMCs. We have also performed multiomic integration to predict microvascular VSMC gene regulatory networks and partitioned SNP heritability across vascular cell types. 

Upon publication, the processed COMMA datasets will be available for interactive exploration through CELLxGENE.

---

## 👩‍💻 Code Structure
<details>
<summary><strong>Click to view full code structure</strong></summary>

**[000_preprocessing](./000_preprocessing/)**

- Code to download and preprocess sc/snRNA-seq and scATAC-seq data. Produces Supplemental Figures 1 & 10.

**[001_scsnRNA-seq](./001_scsnRNA-seq/)**

- Code to benchmark cardiac sc/snRNA-seq integration, re-integrate microvascular cells and annotate cell types. Produces Figure 1, Supplemental Figures 2-5, Supplemental Tables 4-5.

**[002_VSMC_phenotypes](./002_VSMC_phenotypes/)**

- Code to integrate VSMC sc/snRNA-seq data and perform cluster stability analysis, phenotype annotation and comparison with PDGF-BB gene signatures. Produces Figure 2, Supplemental Figure 6, Supplemental Tables 6-7.

**[003_Epicardial_Comparison](./003_Epicardial_Comparison/)**

- Code to integrate COMMA with epicardial coronary artery scRNA-seq data, comparing gene expression and cell type proportions. Produces Figure 3, Supplemental Figure 7, Supplemental Table 9.

**[004_CVD_DGE](./004_CVD_DGE/)**

- Code to perform consensus differential gene expression analysis across cardiovascular disease phenotypes in VSMCs. Produces Figure 4, Supplemental Figures 8-9, Table 1.

**[005_scATAC-seq](./005_scATAC-seq/)**

- Code to integrate cardiac and vascular scATAC-seq data, perform TF prediction, ABC analysis and multiomic integration. Produces Figure 5, Supplemental Figures 11-13, Supplemental Tables 10-11.

**[006_GWAS_Integration](./006_GWAS_Integration/)**

- Code to integrate vascular cell type marker OCRs with known GWAS loci, partitioning CAD heritability and investigating CMD-associated loci. Produces Figure 6, Supplemental Figure 14.

**[CELLxGENE](./CELLxGENE/)**

- Code to format count matrices and fragment files for the CELLxGENE platform.

**[dependencies](./dependencies/)**

- Package versions and dependencies required for different environments.

</details>

---

## 🔁 How to Reproduce

1. Clone the repository:
```bash
git clone https://github.com/tobyb1212/COMMA.git
cd COMMA
```
2. Install required Conda environments, R packages and software [here](./dependencies/).

3. Update absolute file path locations to your directories where relevant.

4. Run analyses in numerical directory and script order where possible.

---

 ##  📋 Metadata

Final metadata for all the sc/snRNA-seq and scATAC-seq data used in this project, will be available upon publication, or from the original studies.

---

 ## 📚 References

The following studies were used in the generation of COMMA. Data from references 14 and 15 were excluded during the generation of COMMA.

1. Litviňuková, M., et al. (2020). Cells of the adult human heart. *Nature*, **588**, 466–472. https://doi.org/10.1038/s41586-020-2797-4
2. Tucker, N. R., et al. (2020). Transcriptional and cellular diversity of the human heart. *Circulation*, **142**, 466–482. https://doi.org/10.1161/CIRCULATIONAHA.119.045401
3. Hocker, J. D., et al. (2021). Cardiac cell type–specific gene regulatory programs and disease risk association. *Science Advances*, **7**, eabf1444. https://doi.org/10.1126/sciadv.abf1444
4. Rao, M., et al. (2021). Resolving the intertwining of inflammation and fibrosis in human heart failure at single-cell level. *Basic Research in Cardiology*, **116**, 55. https://doi.org/10.1007/s00395-021-00897-1
5. Chaffin, M., et al. (2022). Single-nucleus profiling of human dilated and hypertrophic cardiomyopathy. *Nature*, **608**, 174–180. https://doi.org/10.1038/s41586-022-04817-8
6. Koenig, A. L., et al. (2022). Single-cell transcriptomics reveals cell-type-specific diversification in human heart failure. *Nature Cardiovascular Research*, **1**, 263–280. https://doi.org/10.1038/s44161-022-00028-6
7. Amrute, J. M., et al. (2023). Defining cardiac functional recovery in end-stage heart failure at single-cell resolution. *Nature Cardiovascular Research*, **2**, 399–416. https://doi.org/10.1038/s44161-023-00260-8
8. Kanemaru, K., et al. (2023). Spatially resolved multiomics of human cardiac niches. *Nature*, **619**, 801–810. https://doi.org/10.1038/s41586-023-06311-1
9. Read, D. F., et al. (2024). Single-cell analysis of chromatin and expression reveals age- and sex-associated alterations in the human heart. *Communications Biology*, **7**, 1052. https://doi.org/10.1038/s42003-024-06582-y
10. Simonson, B., et al. (2023). Single-nucleus RNA sequencing in ischemic cardiomyopathy reveals common transcriptional profile underlying end-stage heart failure. *Cell Reports*, **42**, 112086. https://doi.org/10.1016/j.celrep.2023.112086
11. Amrute, J. M., et al. (2024). Targeting immune–fibroblast cell communication in heart failure. *Nature*, **635**, 423–433. https://doi.org/10.1038/s41586-024-08008-5
12. Hill, M. C., et al. (2024). Large-scale single-nuclei profiling identifies role for ATRNL1 in atrial fibrillation. *Nature Communications*, **15**, 10002. https://doi.org/10.1038/s41467-024-54296-w
13. Linna-Kuosmanen, S., et al. (2024). Transcriptomic and spatial dissection of human ex vivo right atrial tissue reveals proinflammatory microvascular changes in ischemic heart disease. *Cell Reports Medicine*, **5**, 101556. https://doi.org/10.1016/j.xcrm.2024.101556
14. Selewa, A., et al. (2023). Single-cell genomics improves the discovery of risk variants and genes of atrial fibrillation. *Nature Communications*, **14**, 4999. https://doi.org/10.1038/s41467-023-40505-5
15. Leblanc, F. J. A., et al. (2025). Single-nucleus multi-omics implicates androgen receptor signaling in cardiomyocytes and NR4A1 regulation in fibroblasts during atrial fibrillation. *Nature Cardiovascular Research*, **4**, 433–444. https://doi.org/10.1038/s44161-025-00626-0

---

 ## 📝 License 

This repository is licensed under the [MIT License](./LICENSE).

---

 ## 💻 Contact

For more information, please contact [toby.brown1201@gmail.com](mailto:toby.brown1201@gmail.com).
