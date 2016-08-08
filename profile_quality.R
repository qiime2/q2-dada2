###################################################
# This R script takes a fastq or fastq.gz files
# and outputs a pdf and png of the quality profile of
# the sequences.
#
# Rscript profile_quality.R path/to/fastq/filename.fastq path/to/output/dir
####################################################

# Get character vector of input arguments
# For now: Format is two positional arguments
# path/to/fastq/filename.fastq path/to/output/dir
args <- commandArgs(TRUE)

inp.path <- args[[1]]
out.dir <- args[[2]]
errQuit <- function(mesg) {
  message(mesg)
  q(status=1)
}

# Input file is expected to exist.
if(!file.exists(inp.path)) {
  errQuit("Input filename does not exist.")
}

# Output path is to be a directory and present.
if(!dir.exists(out.dir)) {
  errQuit("Output path does not point to a valid directory.")
}

# Valid input/output -- load libraries and process
library(methods)
library(dada2)
library(ggplot2)

# Plot and save as both pdf and png
p <- plotQualityProfile(inp.path)
ggsave(file.path(out.dir, paste0(basename(inp.path), ".qprofile.pdf")), p,
       device="pdf", width=6, height=4, units="in")
ggsave(file.path(out.dir, paste0(basename(inp.path), ".qprofile.png")), p,
       device="png", width=6, height=4, units="in")
q(status=0)
