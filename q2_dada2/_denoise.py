# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import hashlib

import biom
import skbio
from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import (
    FastqGzFormat, SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)

from q2_dada2._plot import run_commands


def _check_featureless_table(fp):
    with open(fp) as fh:
        # There is a comment line and a header before the feature data
        for line_count, _ in zip(range(3), fh):
            pass
    if line_count < 2:
        raise ValueError("No features remain after denoising. Try adjusting "
                         "your truncation and trim parameter settings.")


def _denoise_helper(cmd, biom_fp, hashed_feature_ids):
    run_commands([cmd])
    _check_featureless_table(biom_fp)
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


def denoise_single(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
                   trunc_len: int, trim_left: int=0, max_ee: float=2.0,
                   trunc_q: int=2, n_threads: int=1,
                   n_reads_learn: int=1000000, hashed_feature_ids: bool=True
                   ) -> (biom.Table, DNAIterator):
    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        cmd = ['run_dada_single.R',
               str(demultiplexed_seqs), biom_fp, temp_dir_name,
               str(trunc_len), str(trim_left), str(max_ee), str(trunc_q),
               str(n_threads), str(n_reads_learn)]
        return _denoise_helper(cmd, biom_fp, hashed_feature_ids)


def denoise_paired(demultiplexed_seqs: SingleLanePerSamplePairedEndFastqDirFmt,
                   trunc_len_f: int, trunc_len_r: int,
                   trim_left_f: int=0, trim_left_r: int=0,
                   max_ee: float=2.0, trunc_q: int=2, n_threads: int=1,
                   n_reads_learn: int=1000000, hashed_feature_ids: bool=True
                   ) -> (biom.Table, DNAIterator):
    with tempfile.TemporaryDirectory() as temp_dir:
        tmp_forward = os.path.join(temp_dir, 'forward')
        tmp_reverse = os.path.join(temp_dir, 'reverse')
        biom_fp = os.path.join(temp_dir, 'output.tsv.biom')
        filt_forward = os.path.join(temp_dir, 'filt_f')
        filt_reverse = os.path.join(temp_dir, 'filt_r')
        for fp in tmp_forward, tmp_reverse, filt_forward, filt_reverse:
            os.mkdir(fp)
        for rp, view in demultiplexed_seqs.sequences.iter_views(FastqGzFormat):
            fp = str(view)
            if 'R1_001.fastq' in rp.name:
                os.link(fp, os.path.join(tmp_forward, rp.name))
            elif 'R2_001.fastq' in rp.name:
                os.link(fp, os.path.join(tmp_reverse, rp.name))

        cmd = ['run_dada_paired.R',
               tmp_forward, tmp_reverse, biom_fp, filt_forward, filt_reverse,
               str(trunc_len_f), str(trunc_len_r),
               str(trim_left_f), str(trim_left_r),
               str(max_ee), str(trunc_q),
               str(n_threads), str(n_reads_learn)]
        return _denoise_helper(cmd, biom_fp, hashed_feature_ids)
