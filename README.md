# R code for the DADA2 QIIME2 plugin

-----------------------------

This repository contains R scripts for use with the [QIIME2 plugin](https://github.com/qiime2/qiime2/wiki/Creating-a-QIIME-2-plugin) version of
DADA2. It is assumed that R and [the dada2 R package](http://benjjneb.github.io/dada2/) have been installed on your
system, and that `Rscript` is in your path.

-----------------------------

`run_dada.R`: This R script takes an input directory of .fastq.gz files
and outputs a tsv file of the dada2 processed sequence
table in "QIIME classic" OTU table format.

`profile_quality.R`: This R script takes a fastq or fastq.gz file 
and outputs a pdf and png of the visualized quality profile of
the sequences.

Example output files are included in the base directory. To replicate those outputs, execute the following shell commands from the base directory of this repository:

```S
Rscript run_dada.R filtered ./seqtab.tsv
Rscript profile_quality.R filtered/F3D0_F_filt.fastq.gz .
```
