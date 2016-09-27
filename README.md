# DADA2 QIIME 2 plugin

This package contains the DADA2 QIIME 2 plugin. It is currently in pre-release status. This documentation is intended for demonstration/testing purposes only. If you're interested in learning to use QIIME 2, you should start at http://2.qiime.org.

## Using DADA2 through QIIME 2

You can install the DADA2 plugin as follows using Miniconda.

```bash
conda create -n q2-dada2 -c r r python=3.5
source activate q2-dada2
conda install -c bioconda bioconductor-dada2
pip install https://github.com/qiime2/qiime2/archive/master.zip https://github.com/qiime2/q2cli/archive/master.zip https://github.com/qiime2/q2-types/archive/master.zip https://github.com/qiime2/q2-feature-table/archive/master.zip https://github.com/benjjneb/q2-dada2/archive/master.zip
```

You can first use DADA2 to explore the quality of a tutorial dataset as follows.

```bash
curl -sO https://dl.dropboxusercontent.com/u/2868868/data/qiime2/artifacts/fmt-tutorial-per-sample-fastq-1p.qza
# this step will take about 2-3 minutes
qiime dada2 plot-qualities --i-demultiplexed-seqs fmt-tutorial-per-sample-fastq-1p.qza --o-visualization quality-plots --p-n 10
qiime tools view quality-plots.qzv
```

You can then use DADA2 to denoise, dereplicate, and remove chimeras from the same dataset as follows.

```bash
# this step will take about 2-3 minutes
qiime dada2 denoise --i-demultiplexed-seqs fmt-tutorial-per-sample-fastq-1p.qza --o-table fmt-tutorial-table.qza --p-trim-left 10 --p-trunc-len 130 --o-representative-sequences fmt-tutorial-rep-seqs.qza
qiime feature-table summarize --i-table fmt-tutorial-table.qza --o-visualization fmt-tutorial-table-summ
qiime tools view fmt-tutorial-table-summ.qzv
```

You can then perform some diversity analyses as follows. First, begin by downloading the sample metadata from [here](https://docs.google.com/spreadsheets/d/16ANHgoFhnpjehCO6ulVPD1b93FDGuDVgA_xh2O4mIRU/edit?usp=sharing) as a ``tsv`` file named ``sample-metadata.tsv``.

```bash
pip install https://github.com/qiime2/q2-feature-table/archive/master.zip https://github.com/qiime2/q2-diversity/archive/master.zip
qiime feature-table rarefy --i-table fmt-tutorial-table.qza --o-rarefied-table fmt-tutorial-table-even200.qza --p-counts-per-sample 200
qiime diversity beta --i-table fmt-tutorial-table-even200.qza --p-metric braycurtis --o-distance-matrix bray-curtis-dm
qiime diversity pcoa --i-distance-matrix bray-curtis-dm.qza --o-pcoa bray-curtis-pc
qiime emperor plot --i-pcoa bray-curtis-pc.qza --m-sample-metadata-file sample-metadata.tsv --o-visualization bray-curtis-emperor
qiime tools view bray-curtis-emperor.qzv
```
