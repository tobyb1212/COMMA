#!/bin/Rscript

library("optparse") #requires this package to parse from the command line
option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL, 
              help="bamfile directory", metavar="character"),
    make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output directory", metavar="character")
); # This is using the directories provided at the command line and turning them into a variable in R. 
#The help and metavar are there to signal what will be output in the help message if needed.
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); # This code basically converts the opt variable into a matrix of values provided at the command line. 
#There are two columns, one for out and one for data.

if (is.null(opt$data)){
  print_help(opt_parser)
  stop("Data directory containing input bamfiles must be supplied", call.=FALSE) #this signifies that the script can't parse the directory from the command line it will return this error message.)
}

if (is.null(opt$out)){
	opt$out <- opt$data
} # If there is no out directory, it will automatically go to the Input Directory

bam.files <- list.files(path = opt$data, pattern = ".bam$", full.names = TRUE) # Defining the string bam.files to the list of files that are in the directory and finish in .bam format.

library(Rsubread) # Loading SubRead from the Bioconductor Bundle

count_fc <- featureCounts(files = bam.files, 
                          annot.ext = '/COMMA/reference_genome/ENSEMBL/GRCr8/Rattus_norvegicus.GRCr8.114.chr.gtf',
                          isGTFAnnotationFile = TRUE,
                          GTF.attrType = 'gene_id',
                          isPairedEnd = TRUE, 
                          strandSpecific = 0, 
                          nthreads = 6,
			  GTF.attrType.extra = c('gene_biotype')			
)

dirname <- opt$out #Setting the Name of the Out Directory 
filename <- paste("rat_pdgf","_count_fc.rds", sep="") #Setting the File Name, Sets the Date as the File Name
saveRDS(count_fc, file = paste(dirname, filename,sep="" )) 
# Saving the Count_fc matrix in the output directory, with the file name followed by suffix, no whitespace between


