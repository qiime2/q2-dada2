# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
from typing import Optional
import hashlib

import biom
import skbio
import qiime2.util
import pandas as pd
import numpy as np

from q2_types.feature_data import DNAIterator, LinkedDNA
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
from q2_dada2._run_dada import _run_dada2


def _check_featureless_table(fp):
    with open(fp) as fh:
        # There is a comment line and a header before the feature data
        for line_count, _ in zip(range(1, 3), fh):
            pass
    if line_count < 2:
        raise ValueError("No features remain after denoising. Try adjusting "
                         "your truncation and trim parameter settings.")


_WHOLE_NUM = (lambda x: x >= 0, 'non-negative')
_NAT_NUM = (lambda x: x > 0, 'greater than zero')
_POOL_STR = (lambda x: x in {'pseudo', 'independent'},
             'pseudo or independent')
_CHIM_STR = (lambda x: x in {'pooled', 'consensus', 'none'},
             'pooled, consensus or none')
_BOOLEAN = (lambda x: type(x) is bool, 'True or False')
# Better to choose to skip, than to implicitly ignore things that KeyError
_SKIP = (lambda x: True, '')
_BOOL = (lambda x: isinstance(x, bool), "Boolean")
_valid_inputs = {
    'trunc_len': _WHOLE_NUM,
    'trunc_len_f': _WHOLE_NUM,
    'trunc_len_r': _WHOLE_NUM,
    'trim_left': _WHOLE_NUM,
    'trim_left_f': _WHOLE_NUM,
    'trim_left_r': _WHOLE_NUM,
    'max_mismatch': _WHOLE_NUM,
    'max_ee': _NAT_NUM,
    'max_ee_f': _NAT_NUM,
    'max_ee_r': _NAT_NUM,
    'trunc_q': _WHOLE_NUM,
    'min_overlap': _WHOLE_NUM,
    'max_merge_mismatch': _WHOLE_NUM,
    'trim_overhang': _BOOLEAN,
    'min_len': _WHOLE_NUM,
    'max_len': _WHOLE_NUM,
    'pooling_method': _POOL_STR,
    'chimera_method': _CHIM_STR,
    'min_fold_parent_over_abundance': _NAT_NUM,
    'allow_one_off': _BOOLEAN,
    'n_threads': _WHOLE_NUM,
    # 0 is technically allowed, but we don't want to support it because it only
    # takes all reads from the first sample (alphabetically by sample id)
    'n_reads_learn': _NAT_NUM,
    # Skipped because they are valid for whole domain of type
    'hashed_feature_ids': _SKIP,
    'demultiplexed_seqs': _SKIP,
    'homopolymer_gap_penalty': _SKIP,
    'band_size': _SKIP,
    'retain_all_samples': _BOOL,
    'retain_unmerged': _BOOL,
    'front': _SKIP,
    'adapter': _SKIP,
    'indels': _SKIP,
}


# TODO: Replace this with Range predicates when interfaces support them better
def _check_inputs(**kwargs):
    for param, arg in kwargs.items():
        check_is_valid, explanation = _valid_inputs[param]
        if not check_is_valid(arg):
            raise ValueError('Argument to %r was %r, should be %s.'
                             % (param, arg, explanation))


def _filepath_to_sample_single(fp):
    return fp.rsplit('_', 4)[0]


def _filepath_to_sample_paired(fp):
    return fp.rsplit('_', 3)[0]


# Since `denoise-single` and `denoise-pyro` are almost identical, break out
# the bulk of the functionality to this helper util. Typechecking is assumed
# to have occurred in the calling functions, this is primarily for making
# sure that DADA2 is able to do what it needs to do.

def _denoise_helper(biom_fp, track_fp, err_track_fp,
                    hashed_feature_ids, retain_all_samples,
                    paired=False, retain_unmerged=False):

    _check_featureless_table(biom_fp)
    with open(biom_fp) as fh:
        table = biom.Table.from_tsv(fh, None, None, None)

    # If we used denoise_paired the barcode was already stripped from the
    # filename to force the files to sort by id and pair up properly
    # see https://github.com/qiime2/q2-dada2/issues/102
    # and https://github.com/qiime2/q2-dada2/pull/125
    filepath_to_sample = _filepath_to_sample_paired if paired \
        else _filepath_to_sample_single

    df = pd.read_csv(track_fp, sep='\t', index_col=0)
    df.index.name = 'sample-id'
    df = df.rename(index=filepath_to_sample)

    PASSED_FILTER = 'percentage of input passed filter'
    NON_CHIMERIC = 'percentage of input non-chimeric'
    CONCATENATED = 'percentage of input concatenated'

    round_cols = {PASSED_FILTER: 2, NON_CHIMERIC: 2}

    df[PASSED_FILTER] = df['filtered'] / df['input'] * 100
    df[NON_CHIMERIC] = df['non-chimeric'] / df['input'] * 100

    col_order = ['input', 'filtered', PASSED_FILTER, 'denoised',
                 'non-chimeric', NON_CHIMERIC]

    # only calculate percentage of input merged if paired end
    if 'merged' in df:
        MERGED = 'percentage of input merged'
        round_cols[MERGED] = 2
        df[MERGED] = df['merged'] / df['input'] * 100
        col_order.insert(4, 'merged')
        col_order.insert(5, MERGED)

    # only calculate percentage of input concatenated if unmerged read pairs
    # were retained
    if 'concatenated' in df:
        round_cols[CONCATENATED] = 2
        df[CONCATENATED] = df['concatenated'] / df['input'] * 100
        insert_at = col_order.index('non-chimeric')
        col_order.insert(insert_at, 'concatenated')
        col_order.insert(insert_at + 1, CONCATENATED)

    # only calculate percentage of input primer-removed if ccs
    if 'primer-removed' in df:
        PASSED_PRIMERREMOVE = 'percentage of input primer-removed'
        round_cols[PASSED_PRIMERREMOVE] = 2
        df[PASSED_PRIMERREMOVE] = df['primer-removed'] / df['input'] * 100
        col_order.insert(1, 'primer-removed')
        col_order.insert(2, PASSED_PRIMERREMOVE)

    df = df[col_order]
    df.fillna(0, inplace=True)
    df = df.round(round_cols)
    metadata = qiime2.Metadata(df)

    # reads in error plot df
    df_err = pd.read_csv(err_track_fp, sep='\t', index_col=0)
    df_err.index.name = 'id'
    df_err.index = df_err.index.astype(str)
    metadata_err = qiime2.Metadata(df_err)

    # Currently the sample IDs in DADA2 are the file names. We make
    # them the sample id part of the filename here.
    sid_map = {id_: filepath_to_sample(id_)
               for id_ in table.ids(axis='sample')}
    table.update_ids(sid_map, axis='sample', inplace=True)
    # Reintroduce empty samples dropped by dada2.
    table_cols = table.ids(axis='observation')
    table_rows = list(set(df.index) - set(table.ids()))

    # We only want to do this if something was actually dropped
    if table_rows:
        table_to_add = biom.Table(np.zeros((len(table_cols), len(table_rows))),
                                  table_cols, table_rows, type="OTU table")
        table = table.concat(table_to_add)
    # This is necessary (instead of just not reintroducing above)
    # dada2 will discard samples which are empty after filtering
    # but will keep samples that are empty after merging
    # so there are potentially samples removed here that were not
    # reintroduced above!
    if not retain_all_samples:
        table = table.remove_empty(axis="sample", inplace=False)

    def _to_sequence(sequence, metadata):
        if retain_unmerged:
            return LinkedDNA(sequence, metadata=metadata)
        return skbio.DNA(sequence, metadata=metadata)

    # The feature IDs in DADA2 are the sequences themselves.
    if hashed_feature_ids:
        # Make feature IDs the md5 sums of the sequences.
        fid_map = {id_: hashlib.md5(id_.encode('utf-8')).hexdigest()
                   for id_ in table.ids(axis='observation')}
        table.update_ids(fid_map, axis='observation', inplace=True)

        rep_sequences = DNAIterator((_to_sequence(k, metadata={'id': v})
                                     for k, v in fid_map.items()))
    else:
        rep_sequences = DNAIterator(
            (_to_sequence(id_, metadata={'id': id_})
             for id_ in table.ids(axis='observation')))

    # initalize and populate DADA2 diagnoistic Stats dictionary
    return table, rep_sequences, metadata, metadata_err


def _denoise_single(demultiplexed_seqs, trunc_len, trim_left, max_ee, trunc_q,
                    max_len, pooling_method, chimera_method,
                    min_fold_parent_over_abundance, allow_one_off,
                    n_threads, n_reads_learn, hashed_feature_ids,
                    homopolymer_gap_penalty, band_size, retain_all_samples):
    _check_inputs(**locals())
    if trunc_len != 0 and trim_left >= trunc_len:
        raise ValueError("trim_left (%r) must be smaller than trunc_len (%r)"
                         % (trim_left, trunc_len))
    if max_len != 0 and max_len < trunc_len:
        raise ValueError("trunc_len (%r) must be no bigger than max_len (%r)"
                         % (trunc_len, max_len))
    # Coerce for single end read analysis
    max_len = 'Inf' if max_len == 0 else max_len

    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        track_fp = os.path.join(temp_dir_name, 'track.tsv')
        err_track_fp = os.path.join(temp_dir_name, 'err_track.tsv')

        _run_dada2(
            input_dir=str(demultiplexed_seqs),
            output_path=str(biom_fp),
            output_track=str(track_fp),
            output_err_track=str(err_track_fp),
            filtered_dir=str(temp_dir_name),
            trunc_len=trunc_len,
            trim_left=trim_left,
            max_ee=max_ee,
            trunc_quality=trunc_q,
            max_len=max_len,
            pooling_method=pooling_method,
            chimera_method=chimera_method,
            min_parental_fold=min_fold_parent_over_abundance,
            allow_one_off=allow_one_off,
            num_threads=n_threads,
            learn_min_reads=n_reads_learn,
            homopolymer_gap_penalty=homopolymer_gap_penalty,
            band_size=band_size
        )

        return _denoise_helper(biom_fp, track_fp, err_track_fp,
                               hashed_feature_ids, retain_all_samples)


def denoise_single(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
                   trunc_len: int, trim_left: int = 0, max_ee: float = 2.0,
                   trunc_q: int = 2, pooling_method: str = 'independent',
                   chimera_method: str = 'consensus',
                   min_fold_parent_over_abundance: float = 1.0,
                   allow_one_off: bool = False,
                   n_threads: int = 1, n_reads_learn: int = 1000000,
                   hashed_feature_ids: bool = True,
                   retain_all_samples: bool = True
                   ) -> (biom.Table, DNAIterator,
                         qiime2.Metadata, qiime2.Metadata):
    return _denoise_single(
        demultiplexed_seqs=demultiplexed_seqs,
        trunc_len=trunc_len,
        trim_left=trim_left,
        max_ee=max_ee,
        trunc_q=trunc_q,
        max_len=0,
        pooling_method=pooling_method,
        chimera_method=chimera_method,
        min_fold_parent_over_abundance=min_fold_parent_over_abundance,
        allow_one_off=allow_one_off,
        n_threads=n_threads,
        n_reads_learn=n_reads_learn,
        hashed_feature_ids=hashed_feature_ids,
        homopolymer_gap_penalty=None,
        band_size=16,
        retain_all_samples=retain_all_samples
    )


def denoise_paired(demultiplexed_seqs: SingleLanePerSamplePairedEndFastqDirFmt,
                   trunc_len_f: int, trunc_len_r: int,
                   trim_left_f: int = 0, trim_left_r: int = 0,
                   max_ee_f: float = 2.0, max_ee_r: float = 2.0,
                   trunc_q: int = 2,
                   min_overlap: int = 12,
                   max_merge_mismatch: int = 0,
                   trim_overhang: bool = False,
                   pooling_method: str = 'independent',
                   chimera_method: str = 'consensus',
                   min_fold_parent_over_abundance: float = 1.0,
                   allow_one_off: bool = False,
                   n_threads: int = 1, n_reads_learn: int = 1000000,
                   hashed_feature_ids: bool = True,
                   retain_all_samples: bool = True,
                   retain_unmerged: bool = False
                   ) -> (biom.Table, DNAIterator,
                         qiime2.Metadata, qiime2.Metadata):
    _check_inputs(**locals())
    if trunc_len_f != 0 and trim_left_f >= trunc_len_f:
        raise ValueError("trim_left_f (%r) must be smaller than trunc_len_f"
                         " (%r)" % (trim_left_f, trunc_len_f))
    if trunc_len_r != 0 and trim_left_r >= trunc_len_r:
        raise ValueError("trim_left_r (%r) must be smaller than trunc_len_r"
                         " (%r)" % (trim_left_r, trunc_len_r))
    with tempfile.TemporaryDirectory() as temp_dir:
        tmp_forward = os.path.join(temp_dir, 'forward')
        tmp_reverse = os.path.join(temp_dir, 'reverse')
        biom_fp = os.path.join(temp_dir, 'output.tsv.biom')
        track_fp = os.path.join(temp_dir, 'track.tsv')
        err_track_fp = os.path.join(temp_dir, 'err_track.tsv')
        filt_forward = os.path.join(temp_dir, 'filt_f')
        filt_reverse = os.path.join(temp_dir, 'filt_r')
        manifest_df = demultiplexed_seqs.manifest.view(pd.DataFrame)

        for fp in tmp_forward, tmp_reverse, filt_forward, filt_reverse:
            os.mkdir(fp)
        for _, fps in manifest_df.iterrows():
            fwd_fp = fps['forward']
            rev_fp = fps['reverse']

            fwd_no_barcode = _remove_barcode(os.path.basename(fps['forward']))
            rev_no_barcode = _remove_barcode(os.path.basename(fps['reverse']))

            qiime2.util.duplicate(fwd_fp, os.path.join(tmp_forward,
                                                       fwd_no_barcode))
            qiime2.util.duplicate(rev_fp, os.path.join(tmp_reverse,
                                                       rev_no_barcode))

        _run_dada2(
            input_dir=str(tmp_forward),
            input_dir_rev=str(tmp_reverse),
            output_path=str(biom_fp),
            output_track=str(track_fp),
            output_err_track=str(err_track_fp),
            filtered_dir=str(filt_forward),
            filtered_dir_rev=str(filt_reverse),
            trunc_len=trunc_len_f,
            trunc_len_rev=trunc_len_r,
            trim_left=trim_left_f,
            trim_left_rev=trim_left_r,
            max_ee=max_ee_f,
            max_ee_rev=max_ee_r,
            trunc_quality=trunc_q,
            min_overlap=min_overlap,
            max_merge_mismatch=max_merge_mismatch,
            pooling_method=pooling_method,
            chimera_method=chimera_method,
            min_parental_fold=min_fold_parent_over_abundance,
            allow_one_off=allow_one_off,
            num_threads=n_threads,
            learn_min_reads=n_reads_learn,
            retain_unmerged=retain_unmerged
        )

        return _denoise_helper(biom_fp, track_fp, err_track_fp,
                               hashed_feature_ids, retain_all_samples,
                               paired=True,
                               retain_unmerged=retain_unmerged)


def _remove_barcode(filename):
    cut = filename.rsplit('_', 3)
    id_ = cut[0].rsplit('_', 1)[0]

    cut = cut[1:]
    cut.insert(0, id_)

    return ('_'.join(cut))


def denoise_pyro(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
                 trunc_len: int, trim_left: int = 0, max_ee: float = 2.0,
                 trunc_q: int = 2, max_len: int = 0,
                 pooling_method: str = 'independent',
                 chimera_method: str = 'consensus',
                 min_fold_parent_over_abundance: float = 1.0,
                 allow_one_off: bool = False,
                 n_threads: int = 1, n_reads_learn: int = 250000,
                 hashed_feature_ids: bool = True,
                 retain_all_samples: bool = True
                 ) -> (biom.Table, DNAIterator,
                       qiime2.Metadata, qiime2.Metadata):
    return _denoise_single(
        demultiplexed_seqs=demultiplexed_seqs,
        trunc_len=trunc_len,
        trim_left=trim_left,
        max_ee=max_ee,
        trunc_q=trunc_q,
        max_len=max_len,
        pooling_method=pooling_method,
        chimera_method=chimera_method,
        min_fold_parent_over_abundance=min_fold_parent_over_abundance,
        allow_one_off=allow_one_off,
        n_threads=n_threads,
        n_reads_learn=n_reads_learn,
        hashed_feature_ids=hashed_feature_ids,
        homopolymer_gap_penalty=1,
        band_size=32,
        retain_all_samples=retain_all_samples)


def denoise_ccs(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
                front: str, adapter: Optional[str] = None,
                max_mismatch: int = 2, indels: bool = False,
                trunc_len: int = 0, trim_left: int = 0, max_ee: float = 2.0,
                trunc_q: int = 2, min_len: int = 20, max_len: int = 0,
                pooling_method: str = 'independent',
                chimera_method: str = 'consensus',
                min_fold_parent_over_abundance: float = 3.5,
                allow_one_off: bool = False,
                n_threads: int = 1, n_reads_learn: int = 1000000,
                hashed_feature_ids: bool = True,
                retain_all_samples: bool = True
                ) -> (biom.Table, DNAIterator,
                      qiime2.Metadata, qiime2.Metadata):
    _check_inputs(**locals())
    if trunc_len != 0 and trim_left >= trunc_len:
        raise ValueError("trim_left (%r) must be smaller than trunc_len (%r)"
                         % (trim_left, trunc_len))
    if max_len != 0 and max_len < trunc_len:
        raise ValueError("trunc_len (%r) must be no bigger than max_len (%r)"
                         % (trunc_len, max_len))
    # Coerce for ccs read analysis
    max_len = 'Inf' if max_len == 0 else max_len

    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        track_fp = os.path.join(temp_dir_name, 'track.tsv')
        err_track_fp = os.path.join(temp_dir_name, 'err_track.tsv')
        nop_fp = os.path.join(temp_dir_name, 'nop')
        filt_fp = os.path.join(temp_dir_name, 'filt')
        for fp in nop_fp, filt_fp:
            os.mkdir(fp)

        _run_dada2(
            input_dir=str(demultiplexed_seqs),
            output_path=str(biom_fp),
            output_track=str(track_fp),
            output_err_track=str(err_track_fp),
            removed_primer_dir=str(nop_fp),
            filtered_dir=str(filt_fp),
            forward_primer=front,
            reverse_primer=adapter,
            max_mismatch=max_mismatch,
            indels=indels,
            trunc_len=trunc_len,
            trim_left=trim_left,
            max_ee=max_ee,
            trunc_quality=trunc_q,
            max_len=max_len,
            min_len=min_len,
            pooling_method=pooling_method,
            chimera_method=chimera_method,
            min_parental_fold=min_fold_parent_over_abundance,
            allow_one_off=allow_one_off,
            num_threads=n_threads,
            learn_min_reads=n_reads_learn,
            band_size=32
        )

        return _denoise_helper(biom_fp, track_fp, err_track_fp,
                               hashed_feature_ids, retain_all_samples)
