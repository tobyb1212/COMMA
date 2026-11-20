library(ArchR)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)

proj <- loadArchRProject(path = '/COMMA/results/scatac/ArchR/initial', showLogo = FALSE)

# Filter Doublets & calculate number of cells filtered
# The main determinant is the filterRatio
proj_doub <- filterDoublets(proj, filterRatio=1.8)
cat(length(getCellNames(ArchRProj = proj))-length(getCellNames(ArchRProj = proj_doub)), "Doublets filtered from the ArchR File")

saveArchRProject(ArchRProj= proj_doub, '/COMMA/results/scatac/ArchR/filtered', load = FALSE)