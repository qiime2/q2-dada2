# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import tempfile
import hashlib

import biom
import skbio
from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt)

from q2_dada2._plot import run_commands


def denoise(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
            trunc_len: int, trim_left: int, max_ee: int=2, truncq: int=2,
            hashed_feature_ids: bool=True) -> (biom.Table, DNAIterator):
    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        cmd = ['run_dada.R', str(demultiplexed_seqs), biom_fp,
               str(trunc_len), str(trim_left), str(max_ee), str(truncq),
               temp_dir_name]
        run_commands([cmd])
        table = biom.Table.from_tsv(open(biom_fp), None, None, None)
        # Currently the sample IDs in DADA2 are the file names. We make
        # them the sample id part of the filename here.
        sid_map = {id_: id_.split('_')[0] for id_ in table.ids(axis='sample')}
        table.update_ids(sid_map, axis='sample', inplace=True)
        # The feature IDs in DADA2 are the sequences themselves.
        if hashed_feature_ids:
            # Make feature IDs the md5 sums of the sequences.
            fid_map = {id_: hashlib.md5(id_.encode('utf-8')).hexdigest()
                       for id_ in table.ids(axis='observation')}
            table.update_ids(fid_map, axis='observation', inplace=True)

            rep_sequences = DNAIterator((skbio.DNA(k, metadata={'id': v})
                                         for k, v in fid_map.items()))
        else:
            rep_sequences = DNAIterator(
                (skbio.DNA(id_, metadata={'id': id_})
                 for id_ in table.ids(axis='observation')))
    return table, rep_sequences
