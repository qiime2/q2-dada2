#!/usr/bin/env Rscript

###################################################
# This R script takes an input two directories of
# .fastq.gz files, corresponding to matched forward
# and reverse sequence files,
# and outputs a tsv file of the dada2 processed sequence
# table. It is intended for use with the QIIME2 plugin
# for DADA2.
#
# Rscript run_dada_faster_paired.R input_dirF input_dirR output.tsv filtered_dirF filtered_dirR 240 160 0 0 2 2 0 100000
####################################################

####################################################
#             DESCRIPTION OF ARGUMENTS             #
####################################################
# NOTE: All numeric arguments should be zero or positive.
# NOTE: All numeric arguments save maxEE are expected to be integers.
# NOTE: Currently the filterered_dirF/R must already exist.
# NOTE: ALL ARGUMENTS ARE POSITIONAL!
#
### FILE SYSTEM ARGUMENTS ###
#
# 1) File path to directory with the FORWARD .fastq.gz files to be processed.
#    Ex: path/to/dir/with/FWD_fastqgzs
#
# 2) File path to directory with the REVERSE .fastq.gz files to be processed.
#    Ex: path/to/dir/with/REV_fastqgzs
#
# 3) File path to output tsv file. If already exists, will be overwritten.
#    Ex: path/to/output_file.tsv
#
# 4) File path to directory to write the filtered FORWARD .fastq.gz files. These files are intermediate
#               for the full workflow. Currently they remain after the script finishes. Directory must
#               already exist.
#    Ex: path/to/dir/with/FWD_fastqgzs/filtered
#
# 5) File path to directory to write the filtered REVERSE .fastq.gz files. These files are intermediate
#               for the full workflow. Currently they remain after the script finishes. Directory must
#               already exist.
#    Ex: path/to/dir/with/REV_fastqgzs/filtered
#
### FILTERING ARGUMENTS ###
#
# 6) truncLenF - The position at which to truncate forward reads. Forward reads shorter
#               than truncLenF will be discarded.
#    Ex: 240
#
# 7) truncLenR - The position at which to truncate reverse reads. Reverse reads shorter
#               than truncLenR will be discarded.
#    Ex: 160
#
# 8) trimLeftF - The number of nucleotides to remove from the start of
#               each forward read. Should be less than truncLenF.
#    Ex: 0
#
# 9) trimLeftR - The number of nucleotides to remove from the start of
#               each reverse read. Should be less than truncLenR.
#    Ex: 0
#
# 10) maxEE - Reads with expected errors higher than maxEE are discarded.
#               Both forward and reverse reads are independently tested.
#    Ex: 2.5
#
# 11) truncQ - Reads are truncated at the first instance of quality score truncQ.
#                If the read is then shorter than truncLen, it is discarded.
#    Ex: 2
#
### SPEED ARGUMENTS ###
#
# 12) nthreads - The number of threads to use.
#                 Special values: 0 - detect available and use all.
#    Ex: 1
#
# 13) nreads_learn - The minimum number of reads to learn the error model from.
#                 Special values: 0 - Use all input reads.
#    Ex: 1000000
#

cat(R.version$version.string, "\n")
args <- commandArgs(TRUE)

inp.dirF <- args[[1]]
inp.dirR <- args[[2]]
out.path <- args[[3]]
filtered.dirF <- args[[4]]
filtered.dirR <- args[[5]]
truncLenF <- as.integer(args[[6]])
truncLenR <- as.integer(args[[7]])
trimLeftF <- as.integer(args[[8]])
trimLeftR <- as.integer(args[[9]])
maxEE <- as.numeric(args[[10]])
truncQ <- as.integer(args[[11]])
nthreads <- as.integer(args[[12]])
nreads.learn <- as.integer(args[[13]])
errQuit <- function(mesg, status=1) {
  message("Error: ", mesg)
  q(status=status)
}

### VALIDATE ARGUMENTS ###

# Input directory is expected to contain .fastq file(s)
# that have not yet been filtered and globally trimmed
# to the same length.
if(!(dir.exists(inp.dirF) && dir.exists(inp.dirR))) {
  errQuit("Input directory does not exist.")
} else {
  unfiltsF <- list.files(inp.dirF, pattern=".fastq.gz$", full.names=TRUE)
  unfiltsR <- list.files(inp.dirR, pattern=".fastq.gz$", full.names=TRUE)
  if(length(unfiltsF) == 0) {
    errQuit("No input forward files with the expected filename format found.")
  }
  if(length(unfiltsR) == 0) {
    errQuit("No input reverse files with the expected filename format found.")
  }
  if(length(unfiltsF) != length(unfiltsR)) {
    errQuit("Different numbers of forward and reverse .fastq.gz files.")
  }
}

# Output path is to be a filename (not a directory) and is to be
# removed and replaced if already present.
if(dir.exists(out.path)) {
  errQuit("Output filename is a directory.")
} else if(file.exists(out.path)) {
  invisible(file.remove(out.path))
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
cat("1) Filtering ")
for(i in seq_along(unfiltsF)) {
  fileNameF = basename(unfiltsF[i])
  filteredFastqF = file.path(filtered.dirF, fileNameF)
  fileNameR = basename(unfiltsR[i])
  filteredFastqR = file.path(filtered.dirR, fileNameR)
  suppressWarnings(fastqPairedFilter(c(unfiltsF[[i]], unfiltsR[[i]]), c(filteredFastqF, filteredFastqR),
                                     truncLen=c(truncLenF, truncLenR), trimLeft=c(trimLeftF, trimLeftR),
                                     maxEE=maxEE, truncQ=truncQ, rm.phix=TRUE))
  if(file.exists(filteredFastqF)) { # Some of the samples reads passed the filter
    cat(".")
  } else {
    cat("x")
  }
}
filtsF <- list.files(filtered.dirF, pattern=".fastq.gz$", full.names=TRUE)
filtsR <- list.files(filtered.dirR, pattern=".fastq.gz$", full.names=TRUE)
cat("\n")
if(length(filtsF) == 0) { # All reads were filtered out
  errQuit("No reads passed the filter (were truncLenF/R longer than the read lengths?)", status=2)
}

### LEARN ERROR RATES ###
# Dereplicate enough samples to get nreads.learn total reads
cat("2) Learning Error Rates\n")
NREADS <- 0
drpsF <- vector("list", length(filtsF))
drpsR <- vector("list", length(filtsR))
for(i in seq_along(filtsF)) {
  drpsF[[i]] <- derepFastq(filtsF[[i]])
  drpsR[[i]] <- derepFastq(filtsR[[i]])
  NREADS <- NREADS + sum(drpsF[[i]]$uniques)
  if(NREADS > nreads.learn) { break }
}
# Run dada in self-consist mode on those samples
drpsF <- drpsF[1:i]
drpsR <- drpsR[1:i]
cat("2a) Forward Reads\n")
ddsF <- dada(drpsF, err=NULL, selfConsist=TRUE, multithread=multithread)
cat("2b) Reverse Reads\n")
ddsR <- dada(drpsR, err=NULL, selfConsist=TRUE, multithread=multithread)
if(i==1) {
  errF <- ddsF$err_out
  errR <- ddsR$err_out
} else {
  errF <- ddsF[[1]]$err_out
}
cat("\n")

### PROCESS ALL SAMPLES ###
# Process samples used to learn error rates
mergers <- vector("list", length(filtsF))
if(i==1) { # breaks list assignment
  mergers[[1]] <- mergePairs(ddsF, drpsF, ddsR, drpsR)
} else {
  mergers[1:i] <- mergePairs(ddsF, drpsF, ddsR, drpsR)
}
rm(drpsF); rm(drpsR); rm(ddsF); rm(ddsR)
# Loop over rest in streaming fashion with learned error rates
cat("3) Denoise remaining samples ")
if(i < length(filtsF)) {
  for(j in seq(i+1,length(filtsF))) {
    drpF <- derepFastq(filtsF[[j]])
    { sink("/dev/null"); ddF <- dada(drpF, err=errF, multithread=multithread); sink(); }
    drpR <- derepFastq(filtsR[[j]])
    { sink("/dev/null"); ddR <- dada(drpR, err=errR, multithread=multithread); sink(); }
    mergers[[j]] <- mergePairs(ddF, drpF, ddR, drpR)
    cat(".")
  }
}
cat("\n")

# Make sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
cat("4) Remove chimeras\n")
seqtab <- removeBimeraDenovo(seqtab, multithread=multithread)

### WRITE OUTPUT AND QUIT ###
# Formatting as tsv plain-text sequence table table
cat("5) Write output\n")
seqtab <- t(seqtab) # QIIME has OTUs as rows
col.names <- basename(filtsF)
col.names[[1]] <- paste0("#OTU ID\t", col.names[[1]])
write.table(seqtab, out.path, sep="\t",
            row.names=TRUE, col.names=col.names, quote=FALSE)
#saveRDS(seqtab, gsub("tsv", "rds", out.path)) ### TESTING
q(status=0)
