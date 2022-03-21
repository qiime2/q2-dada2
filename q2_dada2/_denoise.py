# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import hashlib
import subprocess

import biom
import skbio
import qiime2.util
import pandas as pd

from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import (
    FastqGzFormat, SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)


def run_commands(cmds, verbose=True):
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
    for cmd in cmds:
        if verbose:
            print("\nCommand:", end=' ')
            print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)


def _check_featureless_table(fp):
    with open(fp) as fh:
        # There is a header before the feature data
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
# Better to choose to skip, than to implicitly ignore things that KeyError
_SKIP = (lambda x: True, '')
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
    'min_len': _WHOLE_NUM,
    'max_len': _WHOLE_NUM,
    'pooling_method': _POOL_STR,
    'chimera_method': _CHIM_STR,
    'min_fold_parent_over_abundance': _NAT_NUM,
    'n_threads': _WHOLE_NUM,
    # 0 is technically allowed, but we don't want to support it because it only
    # takes all reads from the first sample (alphabetically by sample id)
    'n_reads_learn': _NAT_NUM,
    # Skipped because they are valid for whole domain of type
    'hashed_feature_ids': _SKIP,
    'demultiplexed_seqs': _SKIP,
    'homopolymer_gap_penalty': _SKIP,
    'band_size': _SKIP,
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


def _filepath_to_sample(fp):
    return fp.rsplit('_', 4)[0]


def _denoise_helper(biom_fp, track_fp, hashed_feature_ids):
    _check_featureless_table(biom_fp)
    with open(biom_fp) as fh:
        table = biom.Table.from_tsv(fh, None, None, None)

    df = pd.read_csv(track_fp, sep='\t', index_col=0)
    df.index.name = 'sample-id'
    df = df.rename(index=_filepath_to_sample)

    PASSED_FILTER = 'percentage of input passed filter'
    NON_CHIMERIC = 'percentage of input non-chimeric'

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

    # Currently the sample IDs in DADA2 are the file names. We make
    # them the sample id part of the filename here.
    sid_map = {id_: _filepath_to_sample(id_)
               for id_ in table.ids(axis='sample')}
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
    return table, rep_sequences, metadata


# Since `denoise-single` and `denoise-pyro` are almost identical, break out
# the bulk of the functionality to this helper util. Typechecking is assumed
# to have occurred in the calling functions, this is primarily for making
# sure that DADA2 is able to do what it needs to do.
def _denoise_single(demultiplexed_seqs, trunc_len, trim_left, max_ee, trunc_q,
                    max_len, pooling_method, chimera_method,
                    min_fold_parent_over_abundance,
                    n_threads, n_reads_learn, hashed_feature_ids,
                    homopolymer_gap_penalty, band_size):
    _check_inputs(**locals())
    if trunc_len != 0 and trim_left >= trunc_len:
        raise ValueError("trim_left (%r) must be smaller than trunc_len (%r)"
                         % (trim_left, trunc_len))
    if max_len != 0 and max_len < trunc_len:
        raise ValueError("trunc_len (%r) must be no bigger than max_len (%r)"
                         % (trunc_len, max_len))
    # Coerce for `run_dada_single.R`
    max_len = 'Inf' if max_len == 0 else max_len

    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        track_fp = os.path.join(temp_dir_name, 'track.tsv')

        cmd = ['run_dada_single.R',
               str(demultiplexed_seqs), biom_fp, track_fp, temp_dir_name,
               str(trunc_len), str(trim_left), str(max_ee), str(trunc_q),
               str(max_len), str(pooling_method), str(chimera_method),
               str(min_fold_parent_over_abundance), str(n_threads),
               str(n_reads_learn), str(homopolymer_gap_penalty),
               str(band_size)]
        try:
            run_commands([cmd])
        except subprocess.CalledProcessError as e:
            if e.returncode == 2:
                raise ValueError(
                    "No reads passed the filter. trunc_len (%r) may be longer"
                    " than read lengths, or other arguments (such as max_ee"
                    " or trunc_q) may be preventing reads from passing the"
                    " filter." % trunc_len)
            else:
                raise Exception("An error was encountered while running DADA2"
                                " in R (return code %d), please inspect stdout"
                                " and stderr to learn more." % e.returncode)
        return _denoise_helper(biom_fp, track_fp, hashed_feature_ids)


def denoise_single(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
                   trunc_len: int, trim_left: int = 0, max_ee: float = 2.0,
                   trunc_q: int = 2, pooling_method: str = 'independent',
                   chimera_method: str = 'consensus',
                   min_fold_parent_over_abundance: float = 1.0,
                   n_threads: int = 1, n_reads_learn: int = 1000000,
                   hashed_feature_ids: bool = True
                   ) -> (biom.Table, DNAIterator, qiime2.Metadata):
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
        n_threads=n_threads,
        n_reads_learn=n_reads_learn,
        hashed_feature_ids=hashed_feature_ids,
        homopolymer_gap_penalty='NULL',
        band_size='16')


def denoise_paired(demultiplexed_seqs: SingleLanePerSamplePairedEndFastqDirFmt,
                   trunc_len_f: int, trunc_len_r: int,
                   trim_left_f: int = 0, trim_left_r: int = 0,
                   max_ee_f: float = 2.0, max_ee_r: float = 2.0,
                   trunc_q: int = 2, min_overlap: int = 12,
                   pooling_method: str = 'independent',
                   chimera_method: str = 'consensus',
                   min_fold_parent_over_abundance: float = 1.0,
                   n_threads: int = 1, n_reads_learn: int = 1000000,
                   hashed_feature_ids: bool = True
                   ) -> (biom.Table, DNAIterator, qiime2.Metadata):
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
        filt_forward = os.path.join(temp_dir, 'filt_f')
        filt_reverse = os.path.join(temp_dir, 'filt_r')
        for fp in tmp_forward, tmp_reverse, filt_forward, filt_reverse:
            os.mkdir(fp)
        for rp, view in demultiplexed_seqs.sequences.iter_views(FastqGzFormat):
            fp = str(view)
            if 'R1_001.fastq' in rp.name:
                qiime2.util.duplicate(fp, os.path.join(tmp_forward, rp.name))
            elif 'R2_001.fastq' in rp.name:
                qiime2.util.duplicate(fp, os.path.join(tmp_reverse, rp.name))

        cmd = ['run_dada_paired.R',
               tmp_forward, tmp_reverse, biom_fp, track_fp, filt_forward,
               filt_reverse,
               str(trunc_len_f), str(trunc_len_r),
               str(trim_left_f), str(trim_left_r),
               str(max_ee_f), str(max_ee_r), str(trunc_q),
               str(min_overlap), str(pooling_method),
               str(chimera_method), str(min_fold_parent_over_abundance),
               str(n_threads), str(n_reads_learn)]
        try:
            run_commands([cmd])
        except subprocess.CalledProcessError as e:
            if e.returncode == 2:
                raise ValueError(
                    "No reads passed the filter. trunc_len_f (%r) or"
                    " trunc_len_r (%r) may be individually longer than"
                    " read lengths, or trunc_len_f + trunc_len_r may be"
                    " shorter than the length of the amplicon + 12"
                    " nucleotides (the length of the overlap). Alternatively,"
                    " other arguments (such as max_ee or trunc_q) may be"
                    " preventing reads from passing the filter."
                    % (trunc_len_f, trunc_len_r))
            else:
                raise Exception("An error was encountered while running DADA2"
                                " in R (return code %d), please inspect stdout"
                                " and stderr to learn more." % e.returncode)
        return _denoise_helper(biom_fp, track_fp, hashed_feature_ids)


def denoise_pyro(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
                 trunc_len: int, trim_left: int = 0, max_ee: float = 2.0,
                 trunc_q: int = 2, max_len: int = 0,
                 pooling_method: str = 'independent',
                 chimera_method: str = 'consensus',
                 min_fold_parent_over_abundance: float = 1.0,
                 n_threads: int = 1, n_reads_learn: int = 250000,
                 hashed_feature_ids: bool = True
                 ) -> (biom.Table, DNAIterator, qiime2.Metadata):
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
        n_threads=n_threads,
        n_reads_learn=n_reads_learn,
        hashed_feature_ids=hashed_feature_ids,
        homopolymer_gap_penalty='-1',
        band_size='32')


def denoise_ccs(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
                front: str, adapter: str,  max_mismatch: int = 2,
                indels: bool = False, trunc_len: int = 0,
                trim_left: int = 0, max_ee: float = 2.0,
                trunc_q: int = 2, min_len: int = 20, max_len: int = 0,
                pooling_method: str = 'independent',
                chimera_method: str = 'consensus',
                min_fold_parent_over_abundance: float = 3.5,
                n_threads: int = 1, n_reads_learn: int = 1000000,
                hashed_feature_ids: bool = True
                ) -> (biom.Table, DNAIterator, qiime2.Metadata):
    _check_inputs(**locals())
    if trunc_len != 0 and trim_left >= trunc_len:
        raise ValueError("trim_left (%r) must be smaller than trunc_len (%r)"
                         % (trim_left, trunc_len))
    if max_len != 0 and max_len < trunc_len:
        raise ValueError("trunc_len (%r) must be no bigger than max_len (%r)"
                         % (trunc_len, max_len))
    # Coerce for `run_dada_ccs.R`
    max_len = 'Inf' if max_len == 0 else max_len

    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        track_fp = os.path.join(temp_dir_name, 'track.tsv')
        nop_fp = os.path.join(temp_dir_name, 'nop')
        filt_fp = os.path.join(temp_dir_name, 'filt')
        for fp in nop_fp, filt_fp:
            os.mkdir(fp)

        cmd = ['run_dada_ccs.R',
               str(demultiplexed_seqs), biom_fp, track_fp, nop_fp, filt_fp,
               str(front), str(adapter), str(max_mismatch), str(indels),
               str(trunc_len), str(trim_left), str(max_ee), str(trunc_q),
               str(min_len), str(max_len), str(pooling_method),
               str(chimera_method), str(min_fold_parent_over_abundance),
               str(n_threads), str(n_reads_learn)]
        try:
            run_commands([cmd])
        except subprocess.CalledProcessError as e:
            if e.returncode == 2:
                raise ValueError(
                    "No reads passed the filter. trunc_len (%r) may be longer"
                    " than read lengths, or other arguments (such as max_ee"
                    " or trunc_q) may be preventing reads from passing the"
                    " filter." % trunc_len)
            else:
                raise Exception("An error was encountered while running DADA2"
                                " in R (return code %d), please inspect stdout"
                                " and stderr to learn more." % e.returncode)
        return _denoise_helper(biom_fp, track_fp, hashed_feature_ids)

def testing_denoise_ccs(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
                front: str, adapter: str,  max_mismatch: int = 2,
                indels: bool = False, trunc_len: int = 0,
                trim_left: int = 0, max_ee: float = 2.0,
                trunc_q: int = 2, min_len: int = 20, max_len: int = 0,
                pooling_method: str = 'independent',
                chimera_method: str = 'consensus',
                min_fold_parent_over_abundance: float = 3.5,
                n_threads: int = 1, n_reads_learn: int = 1000000,
                hashed_feature_ids: bool = True
                ) -> (biom.Table, DNAIterator, qiime2.Metadata):
    _check_inputs(**locals())
    if trunc_len != 0 and trim_left >= trunc_len:
        raise ValueError("trim_left (%r) must be smaller than trunc_len (%r)"
                         % (trim_left, trunc_len))
    if max_len != 0 and max_len < trunc_len:
        raise ValueError("trunc_len (%r) must be no bigger than max_len (%r)"
                         % (trunc_len, max_len))
    # Coerce for `run_dada_ccs.R`
    max_len = 'Inf' if max_len == 0 else max_len

    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        track_fp = os.path.join(temp_dir_name, 'track.tsv')
        nop_fp = os.path.join(temp_dir_name, 'nop')
        filt_fp = os.path.join(temp_dir_name, 'filt')
        for fp in nop_fp, filt_fp:
            os.mkdir(fp)

        cmd = ['run_dada_single_and_ccs.R',
               str(demultiplexed_seqs), biom_fp, track_fp, nop_fp, filt_fp,
               str(front), str(adapter), str(max_mismatch), str(indels), str(trunc_len),
               str(trim_left), str(max_ee), str(trunc_q), str(min_len), str(max_len),
               str(pooling_method), str(chimera_method), str(min_fold_parent_over_abundance), str(n_threads), str(n_reads_learn),
               str('NULL'), str('32')]
        try:
            run_commands([cmd])
        except subprocess.CalledProcessError as e:
            if e.returncode == 2:
                raise ValueError(
                    "No reads passed the filter. trunc_len (%r) may be longer"
                    " than read lengths, or other arguments (such as max_ee"
                    " or trunc_q) may be preventing reads from passing the"
                    " filter." % trunc_len)
            else:
                raise Exception("An error was encountered while running DADA2"
                                " in R (return code %d), please inspect stdout"
                                " and stderr to learn more." % e.returncode)
        return _denoise_helper(biom_fp, track_fp, hashed_feature_ids)

def testing_denoise_single(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
                   trunc_len: int, trim_left: int = 0, max_ee: float = 2.0,
                   trunc_q: int = 2, pooling_method: str = 'independent',
                   chimera_method: str = 'consensus',
                   min_fold_parent_over_abundance: float = 1.0,
                   n_threads: int = 1, n_reads_learn: int = 1000000,
                   hashed_feature_ids: bool = True
                   ) -> (biom.Table, DNAIterator, qiime2.Metadata):
    return testing_secondary_denoise_single(
        demultiplexed_seqs=demultiplexed_seqs,
        trunc_len=trunc_len,
        trim_left=trim_left,
        max_ee=max_ee,
        trunc_q=trunc_q,
        max_len=0,
        pooling_method=pooling_method,
        chimera_method=chimera_method,
        min_fold_parent_over_abundance=min_fold_parent_over_abundance,
        n_threads=n_threads,
        n_reads_learn=n_reads_learn,
        hashed_feature_ids=hashed_feature_ids,
        homopolymer_gap_penalty='NULL',
        band_size='16')


def testing_secondary_denoise_single(demultiplexed_seqs, trunc_len, trim_left, max_ee, trunc_q,
                    max_len, pooling_method, chimera_method,
                    min_fold_parent_over_abundance,
                    n_threads, n_reads_learn, hashed_feature_ids,
                    homopolymer_gap_penalty, band_size):
    _check_inputs(**locals())
    if trunc_len != 0 and trim_left >= trunc_len:
        raise ValueError("trim_left (%r) must be smaller than trunc_len (%r)"
                         % (trim_left, trunc_len))
    if max_len != 0 and max_len < trunc_len:
        raise ValueError("trunc_len (%r) must be no bigger than max_len (%r)"
                         % (trunc_len, max_len))
    # Coerce for `run_dada_single.R`
    max_len = 'Inf' if max_len == 0 else max_len

    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        track_fp = os.path.join(temp_dir_name, 'track.tsv')

        cmd = ['run_dada_single_and_ccs.R',
               str(demultiplexed_seqs), biom_fp, track_fp, str('NULL'),temp_dir_name,
               str('NULL'), str('NULL'), str('NULL'), str('NULL'), str(trunc_len),
               str(trim_left), str(max_ee), str(trunc_q), str('NULL'), str(max_len),
               str(pooling_method), str(chimera_method),str(min_fold_parent_over_abundance), str(n_threads), str(n_reads_learn),
               str(homopolymer_gap_penalty), str(band_size)]
        try:
            run_commands([cmd])
        except subprocess.CalledProcessError as e:
            if e.returncode == 2:
                raise ValueError(
                    "No reads passed the filter. trunc_len (%r) may be longer"
                    " than read lengths, or other arguments (such as max_ee"
                    " or trunc_q) may be preventing reads from passing the"
                    " filter." % trunc_len)
            else:
                raise Exception("An error was encountered while running DADA2"
                                " in R (return code %d), please inspect stdout"
                                " and stderr to learn more." % e.returncode)
        return _denoise_helper(biom_fp, track_fp, hashed_feature_ids)
