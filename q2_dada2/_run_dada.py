# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from rpy2.rinterface import NULL
from rpy2.robjects import default_converter, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import DataFrame as RDataFrame
from rpy2.robjects.vectors import IntVector, ListVector, Matrix as RMatrix
from rpy2.robjects.vectors import StrVector

from q2_dada2._dada_stats._error_model import _error_model_to_dataframe
from q2_dada2._r_utils import _robj_to_pandas_df

dada2 = importr('dada2')


@dataclass(frozen=True)
class _ReadPaths:
    '''
    Sample-aware collection of single- or paired-end FASTQ paths.

    Attributes
    ----------
    sample_ids : tuple[str, ...]
        Sample IDs in manifest-view order.
    forward : StrVector
        Forward FASTQ paths, or the sole direction of single-end reads, in
        manifest-view order.
    reverse : StrVector or None
        Reverse FASTQ paths for paired-end reads in manifest-view order.
    '''
    sample_ids: tuple[str, ...]
    forward: StrVector
    reverse: StrVector | None = None

    @classmethod
    def from_manifest(cls, manifest: pd.DataFrame) -> '_ReadPaths':
        '''
        Build read paths from a single-end, paired-end, or reverse-only
        manifest. For reverse-only single-end data, the reverse paths are
        stored as the primary (``forward``) paths.
        '''
        primary = 'forward' if 'forward' in manifest else 'reverse'
        reverse = None
        if primary == 'forward' and 'reverse' in manifest:
            reverse = StrVector(manifest['reverse'].astype(str).tolist())

        return cls(
            sample_ids=tuple(map(str, manifest.index)),
            forward=StrVector(manifest[primary].astype(str).tolist()),
            reverse=reverse
        )

    def __len__(self) -> int:
        return len(self.sample_ids)

    def output_paths(
        self,
        forward_dir: Path,
        reverse_dir: Path | None = None
    ) -> '_ReadPaths':
        reverse = None
        if self.reverse is not None:
            if reverse_dir is None:
                raise ValueError(
                    'An output directory is required for reverse reads.'
                )
            reverse = StrVector([
                str(reverse_dir / Path(path).name) for path in self.reverse
            ])

        return _ReadPaths(
            sample_ids=self.sample_ids,
            forward=StrVector([
                str(forward_dir / Path(path).name) for path in self.forward
            ]),
            reverse=reverse
        )

    def select(self, keep: Iterable[bool]) -> '_ReadPaths':
        keep = tuple(keep)
        reverse = None
        if self.reverse is not None:
            reverse = StrVector([
                path for path, selected
                in zip(self.reverse, keep, strict=True) if selected
            ])

        return _ReadPaths(
            sample_ids=tuple(
                sample_id for sample_id, selected
                in zip(self.sample_ids, keep, strict=True) if selected
            ),
            forward=StrVector([
                path for path, selected
                in zip(self.forward, keep, strict=True) if selected
            ]),
            reverse=reverse
        )


def _prepare_ccs_reads(
    filtered_dir: Path,
    removed_primer_dir: Path,
    unfiltered: _ReadPaths,
    forward_primer: str,
    reverse_primer: str | None,
    max_mismatch: int,
    indels: bool,
    trunc_len: int,
    trim_left: int,
    max_ee: float,
    trunc_quality: int,
    max_len: int | str,
    min_len: int,
    multithread: bool | int
) -> tuple[_ReadPaths, pd.DataFrame]:
    '''
    Remove primers from CCS reads, then filter and trim the reads.

    Parameters
    ----------
    filtered_dir : Path
        Directory in which to write filtered FASTQ files.
    removed_primer_dir : Path
        Directory in which to write primer-removed FASTQ files.
    unfiltered : _ReadPaths
        Manifest-identified input FASTQ paths and sample IDs.
    forward_primer : str
        Forward primer sequence in the 5' to 3' direction.
    reverse_primer : str or None
        Optional reverse primer sequence in the 5' to 3' direction.
    max_mismatch : int
        Maximum mismatches allowed when matching a primer.
    indels : bool
        Whether insertions and deletions are allowed in primer matches.
    trunc_len : int
        Length at which to truncate reads after primer removal.
    trim_left : int
        Number of bases to remove from the start of each read.
    max_ee : float
        Maximum expected errors allowed in a read.
    trunc_quality : int
        Quality score at which a read is truncated.
    max_len : int or str
        Maximum read length, or ``"Inf"`` for no maximum.
    min_len : int
        Minimum read length after trimming and truncation.
    multithread : bool or int
        Whether to use multiple threads, or the number of threads to use.

    Returns
    -------
    filtered : _ReadPaths
        Filtered FASTQ paths and their sample IDs.
    filtering_stats : pd.DataFrame
        Per-sample input, primer-removed, and filtered read counts.
    '''
    if reverse_primer is None:
        reverse_primer = NULL
    else:
        reverse_primer = dada2.rc(reverse_primer)

    expected_no_primers = unfiltered.output_paths(removed_primer_dir)

    primer_stats_r = dada2.removePrimers(
        fn=unfiltered.forward,
        fout=expected_no_primers.forward,
        primer_fwd=forward_primer,
        primer_rev=reverse_primer,
        max_mismatch=max_mismatch,
        allow_indels=indels,
        orient=True,
        verbose=True
    )

    primer_stats = _robj_to_pandas_df(primer_stats_r)
    primer_stats.index = unfiltered.sample_ids
    no_primers = expected_no_primers.select(primer_stats['reads.out'] > 0)

    if len(no_primers) == 0:
        raise ValueError(
            'No reads passed the Removing Primers step. Did you select the '
            'right primer(s)?'
        )

    expected_filtered = no_primers.output_paths(filtered_dir)

    filtered_stats_r = dada2.filterAndTrim(
        no_primers.forward,
        expected_filtered.forward,
        truncLen=trunc_len,
        trimLeft=trim_left,
        maxEE=max_ee,
        truncQ=trunc_quality,
        rm_phix=False,
        multithread=multithread,
        maxLen=max_len,
        minLen=min_len,
        minQ=3
    )

    filtered_stats = _robj_to_pandas_df(filtered_stats_r)
    filtered_stats.index = no_primers.sample_ids
    filtered = expected_filtered.select(filtered_stats['reads.out'] > 0)
    if len(filtered) == 0:
        raise ValueError(
            'No reads survived filtering and trimming. '
            '(was truncLen longer than the read length?)'
        )

    filtering_stats = primer_stats.join(
        filtered_stats['reads.out'].rename('filtered')
    )

    return filtered, filtering_stats


def _prepare_short_reads(
    filtered_dir: Path,
    unfiltered: _ReadPaths,
    trunc_len: int,
    trim_left: int,
    max_ee: float,
    trunc_quality: int,
    max_len: int | str,
    multithread: bool | int
) -> tuple[_ReadPaths, pd.DataFrame]:
    '''
    Filter and trim single-end or pyrosequencing reads.

    Parameters
    ----------
    filtered_dir : Path
        Directory in which to write filtered forward or single-end FASTQ files.
    unfiltered : _ReadPaths
        Manifest-identified input FASTQ paths and sample IDs.
    trunc_len : int
        Length at which to truncate forward or single-end reads.
    trim_left : int
        Number of bases to remove from the start of each forward or single-end
        read.
    max_ee : float
        Maximum expected errors allowed in a forward or single-end read.
    trunc_quality : int
        Quality score at which reads are truncated.
    max_len : int or str
        Maximum single-end read length, or ``"Inf"`` for no maximum.
    multithread : bool or int
        Whether to use multiple threads, or the number of threads to use.

    Returns
    -------
    filtered : _ReadPaths
        Filtered FASTQ paths and their sample IDs.
    filtering_stats : pd.DataFrame
        Per-sample input and filtered read counts.
    '''
    expected = unfiltered.output_paths(filtered_dir)

    filtering_stats_r = dada2.filterAndTrim(
        unfiltered.forward,
        expected.forward,
        truncLen=trunc_len,
        trimLeft=trim_left,
        maxEE=max_ee,
        truncQ=trunc_quality,
        rm_phix=True,
        multithread=multithread,
        maxLen=max_len
    )

    filtering_stats = _robj_to_pandas_df(filtering_stats_r)
    filtering_stats.index = unfiltered.sample_ids
    filtered = expected.select(filtering_stats['reads.out'] > 0)
    if len(filtered) == 0:
        raise ValueError(
            'No reads survived filtering and trimming. '
            '(was truncLen longer than the read length?)'
        )

    return filtered, filtering_stats


def _prepare_paired_reads(
    filtered_dir: Path,
    filtered_dir_rev: Path,
    unfiltered: _ReadPaths,
    trunc_len: int,
    trunc_len_rev: int,
    trim_left: int,
    trim_left_rev: int,
    max_ee: float,
    max_ee_rev: float,
    trunc_quality: int,
    multithread: bool | int
) -> tuple[_ReadPaths, pd.DataFrame]:
    '''
    Filter paired reads while retaining their manifest sample identities.

    Parameters
    ----------
    filtered_dir : Path
        Directory in which to write filtered forward FASTQ files.
    filtered_dir_rev : Path
        Directory in which to write filtered reverse FASTQ files.
    unfiltered : _ReadPaths
        Manifest-identified forward and reverse input paths.
    trunc_len : int
        Length at which to truncate forward reads.
    trunc_len_rev : int
        Length at which to truncate reverse reads.
    trim_left : int
        Number of bases to remove from the start of each forward read.
    trim_left_rev : int
        Number of bases to remove from the start of each reverse read.
    max_ee : float
        Maximum expected errors allowed in a forward read.
    max_ee_rev : float
        Maximum expected errors allowed in a reverse read.
    trunc_quality : int
        Quality score at which reads are truncated.
    multithread : bool or int
        Whether to use multiple threads, or the number of threads to use.

    Returns
    -------
    filtered : _ReadPaths
        Sample-aware filtered paths for samples with reads passing filtering.
    filtering_stats : pd.DataFrame
        Input and filtered read counts indexed by manifest sample ID.
    '''
    expected = unfiltered.output_paths(filtered_dir, filtered_dir_rev)
    filtering_stats_r = dada2.filterAndTrim(
        unfiltered.forward,
        expected.forward,
        unfiltered.reverse,
        expected.reverse,
        truncLen=[trunc_len, trunc_len_rev],
        trimLeft=[trim_left, trim_left_rev],
        maxEE=[max_ee, max_ee_rev],
        truncQ=trunc_quality,
        rm_phix=True,
        multithread=multithread
    )
    filtering_stats = _robj_to_pandas_df(filtering_stats_r)
    filtering_stats.index = unfiltered.sample_ids

    filtered = expected.select(filtering_stats['reads.out'] > 0)
    if len(filtered) == 0:
        raise ValueError(
            'No reads survived filtering and trimming. '
            '(was truncLen longer than the read length?)'
        )

    return filtered, filtering_stats


@dataclass(frozen=True)
class _ErrorLearningResults:
    '''
    Results from learning one or two (if paired-end reads) DADA2 error models.

    Attributes
    ----------
    forward : ListVector
        Named R list containing `$err_out`, the learned error-rate matrix;
        `$err_in`, the initial error rates; and `$trans`, the observed
        transition counts by nucleotide substitution and quality score.
        See `dada2::learnErrors` for more information.
    reverse : ListVector or None
        Same structure as `forward`, but for reverse reads.
    stats : pd.DataFrame
        Error-model statistics formatted for plotting.
    '''
    forward: ListVector
    reverse: ListVector | None
    stats: pd.DataFrame


@dataclass(frozen=True)
class _Dada2Results:
    '''
    Final in-memory results from a DADA2 workflow.

    Attributes
    ----------
    sequence_table : pd.DataFrame
        Per-sample feature table with samples as rows and sequences as columns.
    filtering_stats : pd.DataFrame
        Per-sample read counts for each completed processing stage.
    error_stats : pd.DataFrame
        Error-model statistics formatted for plotting.
    '''
    sequence_table: pd.DataFrame
    filtering_stats: pd.DataFrame
    error_stats: pd.DataFrame


def _learn_error_models(
    filts: StrVector,
    learn_min_reads: int,
    multithread: bool | int,
    filts_rev: StrVector | None = None,
    pacbio: bool = False,
    homopolymer_gap_penalty: int | None = None,
    band_size: int | None = None
) -> _ErrorLearningResults:
    '''
    Learn DADA2 error models and construct their plotting statistics.

    Parameters
    ----------
    filts : StrVector
        Paths to filtered forward or single-end FASTQ files.
    learn_min_reads : int
        Minimum number of reads used to learn each error model.
    multithread : bool or int
        Whether to use multiple threads, or the number of threads to use.
    filts_rev : StrVector or None
        Paths to filtered reverse FASTQ files for paired-end reads if provided.
    pacbio : bool
        Whether to learn the forward model with DADA2's PacBio error
        estimation function.
    homopolymer_gap_penalty : int or None
        Homopolymer gap penalty used when learning each error model.
    band_size : int or None
        Band size used when learning each error model.

    Returns
    -------
    _ErrorLearningResults
        Learned forward and optional reverse models, plus their combined
        plotting statistics.
    '''
    if pacbio and filts_rev is not None:
        raise ValueError(
            'PacBio error learning does not accept reverse reads.'
        )

    learn_errors_kwargs = {
        'nreads': learn_min_reads,
        'multithread': multithread
    }
    if pacbio:
        learn_errors_kwargs['errorEstimationFunction'] = dada2.PacBioErrfun
        learn_errors_kwargs['BAND_SIZE'] = band_size
    else:
        if homopolymer_gap_penalty is not None:
            learn_errors_kwargs['HOMOPOLYMER_GAP_PENALTY'] = (
                homopolymer_gap_penalty
            )
        if band_size is not None:
            learn_errors_kwargs['BAND_SIZE'] = band_size

    forward = dada2.learnErrors(filts, **learn_errors_kwargs)
    forward_stats = _error_model_to_dataframe(forward)

    if filts_rev is None:
        reverse = None
        stats = forward_stats
    else:
        reverse = dada2.learnErrors(filts_rev, **learn_errors_kwargs)
        reverse_stats = _error_model_to_dataframe(reverse)
        stats = pd.concat(
            [forward_stats.add_prefix('F_'), reverse_stats.add_prefix('R_')],
            axis=1
        )

    return _ErrorLearningResults(
        forward=forward,
        reverse=reverse,
        stats=stats
    )


@dataclass(frozen=True)
class _DenoiseResults:
    '''
    Results from denoising one read direction.

    Attributes
    ----------
    samples : ListVector
        R list containing one `dada-class` object per sample. Each
        `dada-class` object has the following key slots. The native R list is
        retained so it can be passed directly to other DADA2 functions.
            $denoised:
                Integer vector named by inferred sequence and valued by its
                abundance.
            $map:
                Unnamed integer vector of length `derep-class$uniques`.
                Position i is the ith unique sequence. Value at i is the index
                into `dada-class$denoised` that unique sequence i maps (was
                denoised) to. Or, NA if the unique sequence i was removed
                during denoising.

        For the remaining slots see the `dada-class` type in dada2.
    '''
    samples: ListVector

    @classmethod
    def _normalize_denoise_return(
        cls, samples: ListVector
    ) -> '_DenoiseResults':
        '''
        Normalize DADA2's single- and multi-sample return shapes. If a single
        sample is processed a single `dada-class` is returned instead of a
        single-item R list.
        '''
        if 'dada' in samples.rclass:
            normalized = ListVector.from_length(1)
            normalized[0] = samples
            samples = normalized
        return cls(samples=samples)

    @property
    def read_counts(self) -> list[int]:
        return [sum(dada2.getUniques(sample)) for sample in self.samples]


def _denoise_single_reads(
    filts: StrVector,
    err: ListVector,
    pooling_method: str,
    learn_min_reads: int,
    multithread: bool | int,
    homopolymer_gap_penalty: int | None,
    band_size: int | None
) -> _DenoiseResults:
    '''
    Denoise single-end, pyrosequencing, or CCS reads from FASTQ files.

    Parameters
    ----------
    filts : StrVector
        Paths to filtered FASTQ files.
    err : ListVector
        Learned DADA2 error model.
    pooling_method : str
        Pooling method used during denoising. ``"pseudo"`` delegates
        pseudo-pooling to DADA2.
    learn_min_reads : int
        Number of reads supplied to each DADA2 denoising call.
    multithread : bool or int
        Whether to use multiple threads, or the number of threads to use.
    homopolymer_gap_penalty : int or None
        Homopolymer gap penalty passed to DADA2 when provided.
    band_size : int or None
        Band size passed to DADA2 when provided.

    Returns
    -------
    _DenoiseResults
        Per-sample denoised read objects.
    '''
    kwargs = {
        'nreads': learn_min_reads,
        'multithread': multithread,
        'err': err,
        'pool': 'pseudo' if pooling_method == 'pseudo' else False,
        'verbose': False
    }
    if homopolymer_gap_penalty is not None:
        kwargs['HOMOPOLYMER_GAP_PENALTY'] = homopolymer_gap_penalty
    if band_size is not None:
        kwargs['BAND_SIZE'] = band_size

    dds = dada2.dada(filts, **kwargs)

    return _DenoiseResults._normalize_denoise_return(dds)


@dataclass(frozen=True)
class _RetainedUnmergedResults:
    '''
    Retained-unmerged results for one sample.

    Attributes
    ----------
    mergers : pd.DataFrame
        Accepted and retained pairs in DADA2-compatible form.
    concatenated_count : int
        Number of retained unmerged read pairs.
    id_map : pd.DataFrame
        Mapping from temporary sequences to linked sequences.
    '''
    mergers: pd.DataFrame
    concatenated_count: int
    id_map: pd.DataFrame


@dataclass(frozen=True)
class _PairedMergeResults:
    '''
    Results from merging denoised paired-end reads.

    Attributes
    ----------
    merged_reads : list[RDataFrame]
        Per-sample DADA2 merge results used to construct the sequence table.
    merged_counts : list[int]
        Number of merged reads per sample.
    concatenated_counts : list[int]
        Number of retained unmerged reads per sample.
    unmerged_id_map : pd.DataFrame
        Mapping from temporary sequences to linked sequences.
    '''
    merged_reads: list[RDataFrame]
    merged_counts: list[int]
    concatenated_counts: list[int]
    unmerged_id_map: pd.DataFrame


def _denoise_paired_reads(
    reads: _ReadPaths,
    err: ListVector,
    err_rev: ListVector,
    pooling_method: str,
    multithread: bool | int
) -> tuple[_DenoiseResults, _DenoiseResults]:
    '''
    Denoise paired-end reads from FASTQ files.

    Parameters
    ----------
    reads : _ReadPaths
        Sample-aware filtered forward and reverse FASTQ paths.
    err : ListVector
        Learned forward-read DADA2 error model.
    err_rev : ListVector
        Learned reverse-read DADA2 error model.
    pooling_method : str
        Pooling method used during denoising. `"pseudo"` delegates
        pseudo-pooling to DADA2.
    multithread : bool or int
        Whether to use multiple threads, or the number of threads to use.

    Returns
    -------
    denoised_fwd : _DenoiseResults
        Per-sample denoised forward-read objects.
    denoised_rev : _DenoiseResults
        Per-sample denoised reverse-read objects.
    '''
    dds_fwd = dada2.dada(
        reads.forward,
        err=err,
        pool='pseudo' if pooling_method == 'pseudo' else False,
        multithread=multithread,
        verbose=False
    )
    dds_rev = dada2.dada(
        reads.reverse,
        err=err_rev,
        pool='pseudo' if pooling_method == 'pseudo' else False,
        multithread=multithread,
        verbose=False
    )

    return (
        _DenoiseResults._normalize_denoise_return(dds_fwd),
        _DenoiseResults._normalize_denoise_return(dds_rev)
    )


def _retain_unmerged_pairs(
    merged_pairs: pd.DataFrame,
    denoised_fwd: ListVector,
    denoised_rev: ListVector
) -> _RetainedUnmergedResults:
    '''
    Add rejected read pairs to one sample's accepted mergers.

    Parameters
    ----------
    merged_pairs : pd.DataFrame
        DADA2 merge results including accepted and rejected read pairs.
    denoised_fwd : ListVector
        Denoised forward-read object for the sample.
    denoised_rev : ListVector
        Denoised reverse-read object for the sample.

    Returns
    -------
    results : _RetainedUnmergedResults
        Accepted and retained pairs in DADA2-compatible form, the number of
        retained pairs, and their temporary-to-linked sequence mapping.
    '''
    mergers = merged_pairs.loc[
        merged_pairs['accept'], ['sequence', 'abundance']
    ].copy()
    rejected = merged_pairs.loc[
        ~merged_pairs['accept'], ['forward', 'reverse', 'abundance']
    ]
    concatenated_count = int(rejected['abundance'].sum())

    if len(rejected) == 0:
        id_map = pd.DataFrame(
            {
                'temporary': pd.Series(dtype=str),
                'linked': pd.Series(dtype=str)
            }
        )

        return _RetainedUnmergedResults(
            mergers=mergers,
            concatenated_count=concatenated_count,
            id_map=id_map
        )

    sequence_fwd = denoised_fwd.rx2('clustering').rx2('sequence')
    unmerged_fwd = sequence_fwd.rx(
        IntVector(rejected['forward'].to_list())
    )
    sequence_rev = denoised_rev.rx2('clustering').rx2('sequence')
    unmerged_rev = dada2.rc(sequence_rev.rx(
        IntVector(rejected['reverse'].to_list())
    ))

    temporary_sequences = [
        f'{fwd}NNNNNNNNNN{rev}'
        for fwd, rev in zip(unmerged_fwd, unmerged_rev)
    ]
    linked_sequences = [
        f'{fwd} {rev}'
        for fwd, rev in zip(unmerged_fwd, unmerged_rev)
    ]

    id_map = pd.DataFrame(
        {
            'temporary': temporary_sequences,
            'linked': linked_sequences
        }
    )
    retained = pd.DataFrame(
        {
            'sequence': temporary_sequences,
            'abundance': rejected['abundance'].to_list()
        }
    )
    mergers = pd.concat([mergers, retained], ignore_index=True)

    return _RetainedUnmergedResults(
        mergers=mergers,
        concatenated_count=concatenated_count,
        id_map=id_map
    )


def _merge_paired_reads(
    denoised_fwd: _DenoiseResults,
    denoised_rev: _DenoiseResults,
    reads: _ReadPaths,
    min_overlap: int | None,
    max_merge_mismatch: int | None,
    trim_overhang: bool | None,
    retain_unmerged: bool | None
) -> _PairedMergeResults:
    '''
    Merge denoised paired-end reads and optionally retain rejected pairs.

    Parameters
    ----------
    denoised_fwd : _DenoiseResults
        Per-sample denoised forward-read objects.
    denoised_rev : _DenoiseResults
        Per-sample denoised reverse-read objects.
    reads : _ReadPaths
        Sample-aware filtered forward and reverse FASTQ paths.
    min_overlap : int or None
        Minimum overlap required to merge a forward and reverse read.
    max_merge_mismatch : int or None
        Maximum mismatches allowed in the overlap region.
    trim_overhang : bool or None
        Whether to trim overhanging sequence after alignment.
    retain_unmerged : bool or None
        Whether to retain rejected pairs as linked forward and reverse reads.

    Returns
    -------
    results : _PairedMergeResults
        Per-sample merge results, merged and retained-unmerged counts, and the
        temporary-to-linked sequence mapping.
    '''
    merged_counts = []
    concatenated_counts = []
    mergers = []
    mergers_r = []
    unmerged_id_maps = []

    kwargs = {}
    if min_overlap is not None:
        kwargs['minOverlap'] = min_overlap
    if max_merge_mismatch is not None:
        kwargs['maxMismatch'] = max_merge_mismatch
    if trim_overhang is not None:
        kwargs['trimOverhang'] = trim_overhang
    if retain_unmerged is not None:
        kwargs['returnRejects'] = retain_unmerged

    merged_reads_r = dada2.mergePairs(
        denoised_fwd.samples,
        reads.forward,
        denoised_rev.samples,
        reads.reverse,
        **kwargs
    )
    if isinstance(merged_reads_r, RDataFrame):
        merged_reads_r = [merged_reads_r]
    else:
        merged_reads_r = list(merged_reads_r)

    for mp_r, dd_fwd, dd_rev in zip(
            merged_reads_r,
            denoised_fwd.samples,
            denoised_rev.samples,
            strict=True):
        mp = pandas2ri.rpy2py(mp_r)
        if retain_unmerged:
            merged_counts.append(
                int(mp.loc[mp['accept'], 'abundance'].sum())
            )
            retained = _retain_unmerged_pairs(
                merged_pairs=mp,
                denoised_fwd=dd_fwd,
                denoised_rev=dd_rev
            )
            concatenated_counts.append(retained.concatenated_count)
            mergers.append(retained.mergers)
            unmerged_id_maps.append(retained.id_map)
        else:
            mergers_r.append(mp_r)
            merged_counts.append(int(mp['abundance'].sum()))

    if retain_unmerged:
        with localconverter(default_converter + pandas2ri.converter):
            mergers_r = [
                pandas2ri.py2rpy(merger)
                for merger in mergers
            ]

    if unmerged_id_maps:
        unmerged_id_map = pd.concat(unmerged_id_maps, ignore_index=True)
    else:
        unmerged_id_map = pd.DataFrame(
            {
                'temporary': pd.Series(dtype=str),
                'linked': pd.Series(dtype=str)
            }
        )

    return _PairedMergeResults(
        merged_reads=mergers_r,
        merged_counts=merged_counts,
        concatenated_counts=concatenated_counts,
        unmerged_id_map=unmerged_id_map
    )


def _construct_sequence_table(
    samples: ListVector | list[RDataFrame]
) -> RMatrix:
    '''
    Construct a sequence table from per-sample DADA2 results.

    Parameters
    ----------
    samples : ListVector or list[RDataFrame]
        R list of per-sample DADA2 `dada-class` objects for single-end reads,
        or a list of merge-result data frames for paired-end reads.

    Returns
    -------
    RMatrix
        Per-sample sequence table.
    '''
    return dada2.makeSequenceTable(samples)


def _resolve_multithread(num_threads: int | None) -> bool | int:
    '''Convert the requested thread count to DADA2's multithread argument.'''
    if num_threads is None:
        return False
    if num_threads < 0:
        raise ValueError('Number of threads must be a positive number.')
    if num_threads == 0:
        return True
    return num_threads


def _construct_stats_table(
    filtering_stats: pd.DataFrame,
    denoised_counts: Iterable[int],
    non_chimeric_counts: Iterable[int],
    primer_removed: bool = False,
    merged_counts: Iterable[int] | None = None,
    concatenated_counts: Iterable[int] | None = None
) -> pd.DataFrame:
    '''
    Construct the per-sample read tracking table.

    Parameters
    ----------
    filtering_stats : pd.DataFrame
        Per-sample read counts before and after filtering, including counts
        after primer removal when applicable.
    denoised_counts : Iterable[int]
        Number of denoised reads for each sample that passed filtering.
    non_chimeric_counts : Iterable[int]
        Number of non-chimeric reads for each sample that passed filtering.
    primer_removed : bool
        Whether ``filtering_stats`` includes primer-removed read counts.
    merged_counts : Iterable[int] or None
        Number of merged reads for each paired-end sample that passed
        filtering. If omitted, paired-end columns are not added.
    concatenated_counts : Iterable[int] or None
        Number of retained unmerged reads for each paired-end sample that
        passed filtering. If omitted, the concatenated column is not added.

    Returns
    -------
    track : pd.DataFrame
        Per-sample read counts at each applicable processing stage.
    '''
    track = filtering_stats.copy()
    if primer_removed:
        track.columns = ['input', 'primer-removed', 'filtered']
    else:
        track.columns = ['input', 'filtered']

    track['denoised'] = 0
    if merged_counts is not None:
        track['merged'] = 0
    if concatenated_counts is not None:
        track['concatenated'] = 0
    track['non-chimeric'] = 0

    passed_filtering = track['filtered'] > 0
    track.loc[passed_filtering, 'denoised'] = list(denoised_counts)
    if merged_counts is not None:
        track.loc[passed_filtering, 'merged'] = list(merged_counts)
    if concatenated_counts is not None:
        track.loc[passed_filtering, 'concatenated'] = list(
            concatenated_counts
        )
    track.loc[passed_filtering, 'non-chimeric'] = list(non_chimeric_counts)

    return track


def _remove_chimeras(
    sequence_table: RMatrix,
    chimera_method: str,
    min_parental_fold: float,
    allow_one_off: bool,
    multithread: bool | int
) -> RMatrix:
    '''
    Remove chimeric sequences from a DADA2 sequence table.

    Parameters
    ----------
    sequence_table : RMatrix
        DADA2 sequence table to inspect for chimeras.
    chimera_method : str
        Chimera-removal method, or ``"none"`` to skip removal.
    min_parental_fold : float
        Minimum parental abundance fold difference used to identify chimeras.
    allow_one_off : bool
        Whether to identify one-off bimeras as chimeric.
    multithread : bool or int
        Whether to use multiple threads, or the number of threads to use.

    Returns
    -------
    non_chimeric_table : RMatrix
        R sequence table after chimera removal, or the unchanged sequence
        table when chimera removal is skipped.
    '''
    if (
        chimera_method not in {'pooled', 'consensus'}
        or sequence_table.ncol == 0
    ):
        return sequence_table

    return dada2.removeBimeraDenovo(
        sequence_table,
        method=chimera_method,
        minFoldParentOverAbundance=min_parental_fold,
        allowOneOff=allow_one_off,
        multithread=multithread
    )


def _restore_linked_sequences(
    sequence_table: pd.DataFrame,
    unmerged_id_map: pd.DataFrame
) -> pd.DataFrame:
    '''
    Restore retained read pairs to their linked-sequence representation.

    Parameters
    ----------
    sequence_table : pd.DataFrame
        Chimera-filtered sequence table containing temporary sequence IDs.
    unmerged_id_map : pd.DataFrame
        Temporary-to-linked sequence mapping produced while retaining
        unmerged read pairs.

    Returns
    -------
    restored_table : pd.DataFrame
        Sequence table whose retained-pair columns use linked sequences.
    '''
    if len(unmerged_id_map) == 0:
        return sequence_table

    unmerged_id_map = unmerged_id_map.drop_duplicates()
    ambiguous_ids = unmerged_id_map.loc[
        unmerged_id_map['temporary'].duplicated(keep=False), 'temporary'
    ]
    if len(ambiguous_ids) > 0:
        raise ValueError(
            'Unable to uniquely map retained unmerged sequences from the '
            'temporary DADA2-compatible representation to linked sequences.'
        )

    retained_ids = sequence_table.columns.intersection(
        unmerged_id_map['temporary']
    )
    replacements = unmerged_id_map.loc[
        unmerged_id_map['temporary'].isin(retained_ids)
    ].set_index('temporary')['linked'].to_dict()

    return sequence_table.rename(columns=replacements)


def _finalize_dada2_results(
    sequence_table: RMatrix,
    sample_names: Iterable[str],
    filtering_stats: pd.DataFrame,
    error_stats: pd.DataFrame,
    denoised_counts: Iterable[int],
    chimera_method: str,
    min_parental_fold: float,
    allow_one_off: bool,
    multithread: bool | int,
    primer_removed: bool,
    merged_counts: Iterable[int] | None = None,
    concatenated_counts: Iterable[int] | None = None,
    unmerged_id_map: pd.DataFrame | None = None
) -> _Dada2Results:
    '''
    Apply shared post-denoising steps and assemble in-memory results.

    Parameters
    ----------
    sequence_table : RMatrix
        DADA2 sequence table constructed from denoised or merged samples.
    sample_names : Iterable[str]
        Names to assign to sequence-table rows in processing order.
    filtering_stats : pd.DataFrame
        Per-sample filtering counts.
    error_stats : pd.DataFrame
        Error-model statistics formatted for plotting.
    denoised_counts : Iterable[int]
        Number of denoised reads per sample that passed filtering.
    chimera_method : str
        Chimera-removal method, or ``"none"`` to skip removal.
    min_parental_fold : float
        Minimum parental abundance fold difference used to identify chimeras.
    allow_one_off : bool
        Whether to identify one-off bimeras as chimeric.
    multithread : bool or int
        Whether to use multiple threads, or the number of threads to use.
    primer_removed : bool
        Whether primer-removal counts are present in ``filtering_stats``.
    merged_counts : Iterable[int] or None
        Number of merged reads per sample, if paired reads were processed.
    concatenated_counts : Iterable[int] or None
        Number of retained unmerged reads per sample, if requested.
    unmerged_id_map : pd.DataFrame or None
        Temporary-to-linked sequence mapping for retained unmerged reads.

    Returns
    -------
    _Dada2Results
        Final sequence table, read statistics, and error statistics.
    '''
    sequence_table = _remove_chimeras(
        sequence_table=sequence_table,
        chimera_method=chimera_method,
        min_parental_fold=min_parental_fold,
        allow_one_off=allow_one_off,
        multithread=multithread
    )
    sequence_table = _robj_to_pandas_df(sequence_table)

    if unmerged_id_map is not None:
        sequence_table = _restore_linked_sequences(
            sequence_table=sequence_table,
            unmerged_id_map=unmerged_id_map
        )

    sequence_table.index = list(sample_names)

    track = _construct_stats_table(
        filtering_stats=filtering_stats,
        denoised_counts=denoised_counts,
        non_chimeric_counts=sequence_table.sum(axis='columns').values,
        primer_removed=primer_removed,
        merged_counts=merged_counts,
        concatenated_counts=concatenated_counts
    )

    return _Dada2Results(
        sequence_table=sequence_table,
        filtering_stats=track,
        error_stats=error_stats
    )
