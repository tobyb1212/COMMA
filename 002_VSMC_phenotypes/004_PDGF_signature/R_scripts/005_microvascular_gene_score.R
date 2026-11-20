library(edgeR)
library(tidyverse)
library(biomaRt)

setwd("/COMMA/results/pdgf-bb/DGE/")

fc_matrix <-readRDS(file = "count_matrixrat_pdgf_count_fc.rds") # Reading in .RDS Count Matrix

cols <- gsub("\\.bam$", "", grep("^micro_(control|pdgf)", colnames(fc_matrix$counts), value = TRUE))
cols <- cols[order(grepl("pdgf", cols))]
colnames(fc_matrix$counts) <- sub("\\.bam$", "", colnames(fc_matrix$counts))
fc_matrix$counts <- fc_matrix$counts[, cols]

counts <- fc_matrix$counts # Taking Count matrix
genes <- fc_matrix$annotation$GeneID #Taking Gene Names
DGE <- DGEList(counts = counts, genes = genes) # Creating the DGE list Object
colnames(DGE) <- rownames(DGE$samples) 

Treat <- factor(c(rep("Vehicle", 3), rep("PDGF-BB", 3)), levels = c("Vehicle", "PDGF-BB")) # Detailing the Treatment of Samples
data.frame(Sample=colnames(DGE),Treat) # Creating a data frame detailing the Treatment / Biological Study of each sample
design_matrix <- model.matrix(~Treat) # Making the Experimental Design Matrix
rownames(design_matrix) <- colnames(DGE)

# Filter Low Expressed Genes
keep <- filterByExpr(DGE, group=Treat) # Filtering Genes That Have Limited Expression Across Treatment
DGE <- DGE[keep, , keep.lib.sizes=FALSE] # Changing the Library Sizes to Account for This
DGE$samples$lib.size <- colSums(DGE$counts) # Changing the Library Sizes to the New Sum of the Counts

# Normalisation
DGE <- normLibSizes(DGE) # This Minimises the LogFC between genes using the TMM method
DGE$samples # Checking to See Normalisation Factors & If They Make Sense

DGE <- estimateDisp(DGE, design_matrix, robust = TRUE) # Estimating Dispersion Using a GLM Model Based on the Design Matrix
DGE$common.dispersion # This Shows the Value of the common dispesion (BCV)^2
plotBCV(DGE) # This Shows the Tagwise, Trendwise and common BCVs, relative to logCPM. 

fit <- glmQLFit(DGE, design_matrix, robust=TRUE)
qlf <- glmQLFTest(fit) # Performing the QLF F-test
topTags(qlf) # seeing the Top DEGs

tt <- topTags(qlf, n = nrow(qlf), adjust.method = 'BH') # Generating the Top Tags and Adjusting the P-value Using Benjamini-Hochberg
all_genes <- tt$table  # Performing the QLF F-test
deg <- subset(tt$table, FDR < 0.05) # Subsetting the Table And Ecluding Any with FDR Higher Than 0.0
deg_up <- subset(deg, logFC > 1)

rat_mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host = "https://www.ensembl.org")
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")

human_genes <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene",
                "hsapiens_homolog_associated_gene_name"),filters = "ensembl_gene_id",
                values = deg_up$genes,mart = rat_mart)[,3]
human_genes <- unique(na.omit(human_genes))

write.csv(human_genes, "synthetic_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

deg_down <- subset(deg, logFC < -1)

rat_mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host = "https://www.ensembl.org")
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")

human_genes <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene",
                                    "hsapiens_homolog_associated_gene_name"),filters = "ensembl_gene_id",
                     values = deg_down$genes,mart = rat_mart)[,3]
human_genes <- unique(na.omit(human_genes))

write.csv(human_genes, "contractile_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

