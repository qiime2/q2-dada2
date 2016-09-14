# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import hashlib
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
        # Currently the sample IDs in DADA2 are the file names. We make
        # them the sample id part of the filename here.
        sid_map = {id_: id_.split('_')[0] for id_ in table.ids(axis='sample')}
        table.update_ids(sid_map, axis='sample', inplace=True)
        # Currently the feature IDs in DADA2 are the sequences themselves.
        # We make them md5 sums of the sequences here.
        fid_map = {id_: hashlib.md5(id_.encode('utf-8')).hexdigest()
                   for id_ in table.ids(axis='observation')}
        table.update_ids(fid_map, axis='observation', inplace=True)

        rep_sequences = DNAIterator((skbio.DNA(k, metadata={'id': v})
                                     for k, v in fid_map.items()))
    return table, rep_sequences
