# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
from pathlib import Path
from typing import Optional
import hashlib

import biom
import skbio
import qiime2
import pandas as pd
import numpy as np

from q2_types.feature_data import DNAIterator, LinkedDNA
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
from q2_dada2._run_dada import (
    _Dada2Results,
    _construct_sequence_table,
    _denoise_paired_reads,
    _denoise_single_reads,
    _finalize_dada2_results,
    _learn_error_models,
    _merge_paired_reads,
    _ReadPaths,
    _prepare_ccs_reads,
    _prepare_paired_reads,
    _prepare_short_reads,
    _resolve_multithread,
)


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


def _denoise_helper(results: _Dada2Results, hashed_feature_ids,
                    retain_all_samples, retain_unmerged=False):
    if results.sequence_table.shape[1] == 0:
        raise ValueError(
            'No features remain after denoising. Try adjusting your '
            'truncation and trim parameter settings.'
        )

    table = biom.Table(
        results.sequence_table.T.to_numpy(),
        observation_ids=results.sequence_table.columns,
        sample_ids=results.sequence_table.index
    )

    return _assemble_denoise_outputs(
        table=table,
        read_stats=results.filtering_stats.copy(),
        error_stats=results.error_stats.copy(),
        hashed_feature_ids=hashed_feature_ids,
        retain_all_samples=retain_all_samples,
        retain_unmerged=retain_unmerged
    )


def _assemble_denoise_outputs(table, read_stats, error_stats,
                              hashed_feature_ids, retain_all_samples,
                              retain_unmerged=False):
    df = read_stats
    df.index.name = 'sample-id'

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
    df_err = error_stats
    df_err.index.name = 'id'
    df_err.index = df_err.index.astype(str)
    metadata_err = qiime2.Metadata(df_err)

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


# `denoise-single` and `denoise-pyro` differ only in a few DADA2 options, so
# they share the same composed single-read workflow.
def _denoise_single_or_pyro(
    demultiplexed_seqs, trunc_len, trim_left, max_ee, trunc_q,
    max_len, pooling_method, chimera_method,
    min_fold_parent_over_abundance, allow_one_off,
    n_threads, n_reads_learn, hashed_feature_ids,
    homopolymer_gap_penalty, band_size, retain_all_samples
):
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
        temp_dir = Path(temp_dir_name)
        multithread = _resolve_multithread(n_threads)
        unfiltered = _ReadPaths.from_manifest(
            demultiplexed_seqs.manifest.view(pd.DataFrame)
        )
        filtered, filtering_stats = _prepare_short_reads(
            filtered_dir=temp_dir,
            unfiltered=unfiltered,
            trunc_len=trunc_len,
            trim_left=trim_left,
            max_ee=max_ee,
            trunc_quality=trunc_q,
            max_len=max_len,
            multithread=multithread
        )
        error_models = _learn_error_models(
            filts=filtered.forward,
            filts_rev=None,
            learn_min_reads=n_reads_learn,
            multithread=multithread,
            pacbio=False,
            homopolymer_gap_penalty=homopolymer_gap_penalty,
            band_size=band_size
        )
        denoised = _denoise_single_reads(
            filts=filtered.forward,
            err=error_models.forward,
            pooling_method=pooling_method,
            learn_min_reads=n_reads_learn,
            multithread=multithread,
            homopolymer_gap_penalty=homopolymer_gap_penalty,
            band_size=band_size
        )
        sequence_table = _construct_sequence_table(denoised.samples)
        results = _finalize_dada2_results(
            sequence_table=sequence_table,
            sample_names=filtered.sample_ids,
            filtering_stats=filtering_stats,
            error_stats=error_models.stats,
            denoised_counts=denoised.read_counts,
            chimera_method=chimera_method,
            min_parental_fold=min_fold_parent_over_abundance,
            allow_one_off=allow_one_off,
            multithread=multithread,
            primer_removed=False
        )

        return _denoise_helper(
            results=results,
            hashed_feature_ids=hashed_feature_ids,
            retain_all_samples=retain_all_samples
        )


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
    return _denoise_single_or_pyro(
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
    return _denoise_single_or_pyro(
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
        homopolymer_gap_penalty=-1,
        band_size=32,
        retain_all_samples=retain_all_samples)


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

    with tempfile.TemporaryDirectory() as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        filt_forward = temp_dir / 'filt_f'
        filt_reverse = temp_dir / 'filt_r'
        manifest_df = demultiplexed_seqs.manifest.view(pd.DataFrame)

        for directory in (filt_forward, filt_reverse):
            directory.mkdir()

        unfiltered = _ReadPaths.from_manifest(manifest_df)

        multithread = _resolve_multithread(n_threads)
        filtered, filtering_stats = _prepare_paired_reads(
            filtered_dir=filt_forward,
            filtered_dir_rev=filt_reverse,
            unfiltered=unfiltered,
            trunc_len=trunc_len_f,
            trunc_len_rev=trunc_len_r,
            trim_left=trim_left_f,
            trim_left_rev=trim_left_r,
            max_ee=max_ee_f,
            max_ee_rev=max_ee_r,
            trunc_quality=trunc_q,
            multithread=multithread
        )

        error_models = _learn_error_models(
            filts=filtered.forward,
            filts_rev=filtered.reverse,
            learn_min_reads=n_reads_learn,
            multithread=multithread
        )
        if error_models.reverse is None:
            raise RuntimeError(
                'Paired error learning returned no reverse model.'
            )

        denoised_fwd, denoised_rev = _denoise_paired_reads(
            reads=filtered,
            err=error_models.forward,
            err_rev=error_models.reverse,
            pooling_method=pooling_method,
            multithread=multithread
        )
        merged = _merge_paired_reads(
            denoised_fwd=denoised_fwd,
            denoised_rev=denoised_rev,
            reads=filtered,
            min_overlap=min_overlap,
            max_merge_mismatch=max_merge_mismatch,
            trim_overhang=trim_overhang,
            retain_unmerged=retain_unmerged
        )
        sequence_table = _construct_sequence_table(merged.merged_reads)
        results = _finalize_dada2_results(
            sequence_table=sequence_table,
            sample_names=filtered.sample_ids,
            filtering_stats=filtering_stats,
            error_stats=error_models.stats,
            denoised_counts=denoised_fwd.read_counts,
            chimera_method=chimera_method,
            min_parental_fold=min_fold_parent_over_abundance,
            allow_one_off=allow_one_off,
            multithread=multithread,
            primer_removed=False,
            merged_counts=merged.merged_counts,
            concatenated_counts=(
                merged.concatenated_counts if retain_unmerged else None
            ),
            unmerged_id_map=merged.unmerged_id_map
        )

        return _denoise_helper(
            results=results,
            hashed_feature_ids=hashed_feature_ids,
            retain_all_samples=retain_all_samples,
            retain_unmerged=retain_unmerged
        )


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
        temp_dir = Path(temp_dir_name)
        removed_primer_dir = temp_dir / 'nop'
        filtered_dir = temp_dir / 'filt'
        removed_primer_dir.mkdir()
        filtered_dir.mkdir()

        multithread = _resolve_multithread(n_threads)
        unfiltered = _ReadPaths.from_manifest(
            demultiplexed_seqs.manifest.view(pd.DataFrame)
        )
        filtered, filtering_stats = _prepare_ccs_reads(
            filtered_dir=filtered_dir,
            removed_primer_dir=removed_primer_dir,
            unfiltered=unfiltered,
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
            multithread=multithread
        )
        error_models = _learn_error_models(
            filts=filtered.forward,
            filts_rev=None,
            learn_min_reads=n_reads_learn,
            multithread=multithread,
            pacbio=True,
            band_size=32
        )
        denoised = _denoise_single_reads(
            filts=filtered.forward,
            err=error_models.forward,
            pooling_method=pooling_method,
            learn_min_reads=n_reads_learn,
            multithread=multithread,
            homopolymer_gap_penalty=None,
            band_size=32
        )
        sequence_table = _construct_sequence_table(denoised.samples)
        results = _finalize_dada2_results(
            sequence_table=sequence_table,
            sample_names=filtered.sample_ids,
            filtering_stats=filtering_stats,
            error_stats=error_models.stats,
            denoised_counts=denoised.read_counts,
            chimera_method=chimera_method,
            min_parental_fold=min_fold_parent_over_abundance,
            allow_one_off=allow_one_off,
            multithread=multithread,
            primer_removed=True
        )

        return _denoise_helper(
            results=results,
            hashed_feature_ids=hashed_feature_ids,
            retain_all_samples=retain_all_samples
        )
