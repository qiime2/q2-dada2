#!/usr/bin/env Rscript

###################################################
# This R script takes two input directories of .fastq.gz files,
# which are assumed to be paired forward/reverse read files,
# and outputs a tsv file of the dada2 processed sequence
# table. It is intended for use with the QIIME2 plugin
# for DADA2.
####################################################

# Get character vector of input arguments
args <- commandArgs(TRUE)

inp.dirF <- args[[1]]
inp.dirR <- args[[2]]
out.path <- args[[3]]
truncLenF <- as.integer(args[[4]])
truncLenR <- as.integer(args[[5]])
filteredFastqOutputDirF <- args[[6]]
filteredFastqOutputDirR <- args[[7]]
# BJC: Revisit this. Does QIIME need to control formation of temp directories?
errQuit <- function(mesg) {
  message(mesg)
  q(status=1)
}

# Input directories are expected to contain .fastq file(s)
# that have not yet been filtered and globally trimmed
# to the same length.
if(!dir.exists(inp.dirF) || !dir.exists(inp.dirR)) {
  errQuit("Input directory does not exist.")
} else {
  unfiltsF <- list.files(inp.dirF, pattern=".fastq.gz$", full.names=TRUE)
  if(length(unfiltsF) == 0) {
    errQuit("No input files with the expected filename format found.")
  }
  unfiltsR <- list.files(inp.dirR, pattern=".fastq.gz$", full.names=TRUE)
  if(length(unfiltsR) == 0) {
    errQuit("No input files with the expected filename format found.")
  }
  if(length(unfiltsF) != length(unfiltsR)) {
    errQuit("Different numbers of forward and reverse .fastq.gz files.")
  }
}

# Output path is to be a filename (not a directory) and is to be
# removed and replaced if already present.
if(dir.exists(out.path)) {
  errQuit("Output filename points to pre-existing directory.")
} else if(file.exists(out.path)) {
  file.remove(out.path)
}

# Valid input/output -- load libraries and process
library(methods)
library(dada2)

# Trim and filter
# This is adapted from the example provided in the DADA2 tutorial.
for(i in seq_along(unfiltsF)) {
  fileNameF = basename(unfiltsF[i])
  filteredFastqF = file.path(filteredFastqOutputDirF, fileNameF)
  fileNameR = basename(unfiltsR[i])
  filteredFastqR = file.path(filteredFastqOutputDirR, fileNameR)
  fastqPairedFilter(c(unfiltsF[i], unfiltsR[i]), c(filteredFastqF, filteredFastqR),
                    truncLen=c(truncLenF, truncLenR), rm.phix=TRUE)
}
filtsF <- list.files(filteredFastqOutputDirF, pattern=".fastq.gz$", full.names=TRUE)
filtsR <- list.files(filteredFastqOutputDirR, pattern=".fastq.gz$", full.names=TRUE)

## Forward
# Dereplicate
drpsF <- derepFastq(filtsF)
# Sample inference
ddsF <- dada(drpsF, err=NULL, selfConsist = TRUE)
## Reverse
# Dereplicate
drpsR <- derepFastq(filtsR)
# Sample inference
ddsR <- dada(drpsR, err=NULL, selfConsist=TRUE)

# Merge
mergers <- mergePairs(ddsF, drpsF, ddsR, drpsR)

# Make sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab)

# Write output and quit
# Formatting as tsv plain-text OTU table
seqtab <- t(seqtab) # QIIME has OTUs as rows
col.names <- basename(filts)
col.names[[1]] <- paste0("#OTU ID\t", col.names[[1]])
write.table(seqtab, out.path, sep="\t",
            row.names=TRUE, col.names=col.names, quote=FALSE)
q(status=0)
