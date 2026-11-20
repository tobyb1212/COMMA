library(ArchR)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)

proj <- loadArchRProject(path = '/COMMA/results/scatac/ArchR/initial', showLogo = FALSE)

colnames(proj@cellColData)
colnames(proj@sampleColData)

# Supplemental Table 9
obs <- readxl::read_xlsx('/COMMA/results/scrna/metadata/FINAL_METADATA.xlsx', col_names = TRUE, sheet = 2)
obs <- obs %>% select(-c(...21, ...22))

new_sampleColData <- obs[match(rownames(proj@sampleColData), obs$fastq_name), ]

for (c in colnames(new_sampleColData)){
  proj = addSampleColData(ArchRProj = proj, data = new_sampleColData[[c]], sample = unique(proj$Sample), name = c)
}

new_cellColData <- as.data.frame(proj@cellColData) %>%
    rownames_to_column("cellNames") %>%
    select(cellNames, Sample) %>%
    merge(obs, by.x = "Sample", by.y = "fastq_name", all.x = TRUE) %>%
    select(-Sample) %>%
    column_to_rownames("cellNames")

new_cellColData <- new_cellColData[rownames(proj@cellColData), ]

for (c in colnames(new_cellColData)) {
  proj <- addCellColData(ArchRProj = proj, data = new_cellColData[[c]], name = c, cells = getCellNames(proj))
}

colnames(proj@cellColData)
colnames(proj@sampleColData)

saveArchRProject(ArchRProj = proj,'/COMMA/results/scatac/ArchR/initial', load = TRUE)