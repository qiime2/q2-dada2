#!/usr/bin/env Rscript

###################################################
# This R script takes an input directory of .fastq.gz files
# and outputs a tsv file of the dada2 processed sequence
# table. It is intended for use with the QIIME2 plugin
# for DADA2.
#
# Rscript run_dada2.R input_dir output.tsv truncLen trimLeft maxEE truncQ filtered_dir
####################################################

####################################################
#             DESCRIPTION OF ARGUMENTS             #
####################################################
# NOTE: ALL ARGUMENTS ARE POSITIONAL!
#
# 1) File path to directory with the .fastq.gz files to be processed.
#    Ex: path/to/dir/with/fastqgzs
#
# 2) File path to output tsv file. If already exists, will be overwritten.
#    Ex: path/to/output_file.tsv
#
# 3) truncLen - The position at which to truncate reads. Reads shorter
#               than truncLen will be discarded.
#    Ex: 150
#
# 4) trimLeft - The number of nucleotides to remove from the start of
#               each read. Should be less than truncLen for obvious reasons.
#    Ex: 0
#
# 5) maxEE - Reads with expected errors higher than maxEE are discarded.
#    Ex: 3
#
# 6) truncQ - Reads are truncated at the first instance of quality score truncQ.
#                If the read is then shorter than truncLen, it is discarded.
#    Ex: 2
#
# 7) File path to write the filtered .fastq.gz files. These files remain after
#               the script completes.
#    Ex: path/to/dir/with/fastqgzs/filtered
#


cat(R.version$version.string, "\n")
args <- commandArgs(TRUE)

inp.dir <- args[[1]]
out.path <- args[[2]]
truncLen <- as.integer(args[[3]])
trimLeft <- as.integer(args[[4]])
maxEE <- as.integer(args[[5]])
truncQ <- as.integer(args[[6]])
filteredFastqOutputDir <- args[[7]]
errQuit <- function(mesg) {
  message(mesg)
  q(status=1)
}

# Input directory is expected to contain .fastq file(s)
# that have not yet been filtered and globally trimmed
# to the same length.
if(!dir.exists(inp.dir)) {
  errQuit("Input directory does not exist.")
} else {
  unfilts <- list.files(inp.dir, pattern=".fastq.gz$", full.names=TRUE)
  if(length(unfilts) == 0) {
    errQuit("No input files with the expected filename format found.")
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
cat("DADA2 R package version:", as.character(packageVersion("dada2")), "\n")

# Trim and filter
# This is adapted from the example provided in the DADA2 tutorial.
for(i in seq_along(unfilts)) {
  fileName = basename(unfilts[i])
  filteredFastq = file.path(filteredFastqOutputDir, fileName)
  fastqFilter(unfilts[i], filteredFastq, truncLen=truncLen, trimLeft=trimLeft,
              maxEE=maxEE, truncQ=truncQ, rm.phix=TRUE)
}
filts <- list.files(filteredFastqOutputDir, pattern=".fastq.gz$",
                    full.names=TRUE)

# Dereplicate
drps <- derepFastq(filts)

# Sample inference
dds <- dada(drps, err=NULL, selfConsist = TRUE)

# Make sequence table
seqtab <- makeSequenceTable(dds)

# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab)

# Write output and quit
# Formatting as tsv plain-text OTU table
seqtab <- t(seqtab) # QIIME has OTUs as rows
col.names <- basename(filts)
col.names[[1]] <- paste0("#OTU ID\t", col.names[[1]])
write.table(seqtab, out.path, sep="\t",
            row.names=TRUE, col.names=col.names, quote=FALSE)
saveRDS(seqtab, gsub("tsv", "rds", out.path)) ### TESTING
q(status=0)
