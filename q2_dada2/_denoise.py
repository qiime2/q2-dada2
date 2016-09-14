# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import tempfile
import subprocess

import biom
import skbio
from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt)


def denoise(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
            trunc_len: int, trim_left:int) -> (biom.Table, DNAIterator):
    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        cmd = ['run_dada.R', str(demultiplexed_seqs), biom_fp,
               str(trunc_len), str(trim_left), temp_dir_name]
        subprocess.run(cmd, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        table = biom.Table.from_tsv(open(biom_fp), None, None, None)
        # Currently the feature IDs in DADA2 are the sequences themselves,
        # so this result has the sequence both as the id and the data. We're
        # ultimately going to change the ids to be something else (see
        # https://github.com/benjjneb/q2-dada2/issues/4).
        rep_sequences = DNAIterator((skbio.DNA(id_, metadata={'id': id_})
                                     for id_ in table.ids(axis='observation')))
    return table, rep_sequences
