#!/usr/bin/env Rscript

###################################################
# This R script takes an input directory of .fastq.gz files
# and outputs a tsv file of the dada2 processed sequence
# table. It is intended for use with the QIIME2 plugin
# for DADA2.
#
# Rscript run_dada_faster.R input_dir output.tsv truncLen trimLeft filtered_dir
####################################################

####################################################
#             DESCRIPTION OF ARGUMENTS             #
####################################################
# NOTE: All numeric arguments should be zero or positive.
# NOTE: All numeric arguments save maxEE are expected to be integers.
# NOTE: Currently the filterered_dir must already exist.
# NOTE: ALL ARGUMENTS ARE POSITIONAL!
#
### FILE SYSTEM ARGUMENTS ###
#
# 1) File path to directory with the .fastq.gz files to be processed.
#    Ex: path/to/dir/with/fastqgzs
#
# 2) File path to output tsv file. If already exists, will be overwritten.
#    Ex: path/to/output_file.tsv
#
# 3) File path to write the filtered .fastq.gz files. These files are interemediate
#               for the full workflow. Currently they remain after the script finishes.
#    Ex: path/to/dir/with/fastqgzs/filtered
#
### FILTERING ARGUMENTS ###
# 
# 4) truncLen - The position at which to truncate reads. Reads shorter
#               than truncLen will be discarded.
#    Ex: 150
#
# 5) trimLeft - The number of nucleotides to remove from the start of
#               each read. Should be less than truncLen for obvious reasons.
#    Ex: 0
#
# 6) maxEE - Reads with expected errors higher than maxEE are discarded.
#    Ex: 2.5
#
# 7) truncQ - Reads are truncated at the first instance of quality score truncQ.
#                If the read is then shorter than truncLen, it is discarded.
#    Ex: 2
#
### SPEED ARGUMENTS ###
#
# 8) nthreads - The number of threads to use.
#                 Special values: 0 - detect available and use all.
#    Ex: 1
#
# 9) nreads_learn - The minimum number of reads to learn the error model from.
#                 Special values: 0 - Use all input reads.
#    Ex: 1000000
#

cat(R.version$version.string, "\n")
args <- commandArgs(TRUE)

inp.dir <- args[[1]]
out.path <- args[[2]]
filtered.dir <- args[[3]]
truncLen <- as.integer(args[[4]])
trimLeft <- as.integer(args[[5]])
maxEE <- as.numeric(args[[6]])
truncQ <- as.integer(args[[7]])
nthreads <- as.integer(args[[8]])
nreads.learn <- as.integer(args[[9]])
errQuit <- function(mesg) {
  message(mesg)
  q(status=1)
}

### VALIDATE ARGUMENTS ###

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

# Convert nthreads to the logical/numeric expected by dada2
if(nthreads < 0) { 
  errQuit("nthreads must be non-negative.") 
} else if(nthreads == 0) { 
  multithread <- TRUE # detect and use all
} else if(nthreads == 1) {
  multithread <- FALSE 
} else { 
  multithread <- nthreads
}

### LOAD LIBRARIES ###
library(methods)
library(dada2)
cat("DADA2 R package version:", as.character(packageVersion("dada2")), "\n")

### TRIM AND FILTER ###
# This is adapted from the example provided in the DADA2 tutorial.
for(i in seq_along(unfilts)) {
  fileName = basename(unfilts[i])
  filteredFastq = file.path(filtered.dir, fileName)
  fastqFilter(unfilts[i], filteredFastq, truncLen=truncLen, trimLeft=trimLeft,
              maxEE=maxEE, truncQ=truncQ, rm.phix=TRUE)
}
filts <- list.files(filtered.dir, pattern=".fastq.gz$",
                    full.names=TRUE)

### RUN DADA2 PIPELINE ###
# Dereplicate
drps <- derepFastq(filts)

# Sample inference
dds <- dada(drps, err=NULL, selfConsist = TRUE, multithread=multithread)

# Make sequence table
seqtab <- makeSequenceTable(dds)

# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab)

### WRITE OUTPUT AND QUIT ###
# Formatting as tsv plain-text sequence table table
seqtab <- t(seqtab) # QIIME has OTUs as rows
col.names <- basename(filts)
col.names[[1]] <- paste0("#OTU ID\t", col.names[[1]])
write.table(seqtab, out.path, sep="\t",
            row.names=TRUE, col.names=col.names, quote=FALSE)
q(status=0)
