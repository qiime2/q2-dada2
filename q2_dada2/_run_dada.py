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


def _find_fastq_files(directory: Path) -> list[Path]:
    # mimic R's list.files(pattern=...) behaviour which excludes dotfiles
    return sorted(
        path for path in directory.glob('*.fastq.gz')
        if not path.name.startswith('.')
    )


def _prepare_ccs_reads(
    filtered_dir: Path,
    removed_primer_dir: Path,
    unfilts: list[Path],
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
) -> tuple[StrVector, pd.DataFrame]:
    '''
    Remove primers from CCS reads, then filter and trim the reads.

    Parameters
    ----------
    filtered_dir : Path
        Directory in which to write filtered FASTQ files.
    removed_primer_dir : Path
        Directory in which to write primer-removed FASTQ files.
    unfilts : list[Path]
        Input FASTQ files from which primers will be removed.
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
    filts : StrVector
        Paths to the filtered FASTQ files.
    filtering_stats : pd.DataFrame
        Per-sample input, primer-removed, and filtered read counts.
    '''
    if reverse_primer is None:
        reverse_primer = NULL
    else:
        reverse_primer = dada2.rc(reverse_primer)

    removed_primers = [removed_primer_dir / path.name for path in unfilts]

    no_primers = dada2.removePrimers(
        fn=StrVector([str(path) for path in unfilts]),
        fout=StrVector([str(path) for path in removed_primers]),
        primer_fwd=forward_primer,
        primer_rev=reverse_primer,
        max_mismatch=max_mismatch,
        allow_indels=indels,
        orient=True,
        verbose=True
    )

    removed_primers = _find_fastq_files(removed_primer_dir)

    if len(removed_primers) == 0:
        raise ValueError(
            'No reads passed the Removing Primers step. Did you select the '
            'right primer(s)?'
        )

    filts = StrVector([
        str(filtered_dir / path.name) for path in removed_primers
    ])

    filtered_out = dada2.filterAndTrim(
        StrVector([str(path) for path in removed_primers]),
        filts,
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

    filts = StrVector([
        str(path) for path in _find_fastq_files(filtered_dir)
    ])
    if len(filts) == 0:
        raise ValueError(
            'No reads survived filtering and trimming. '
            '(was truncLen longer than the read length?)'
        )

    filtering_stats = pd.concat(
        [
            _robj_to_pandas_df(no_primers),
            _robj_to_pandas_df(filtered_out)['reads.out']
        ],
        axis=1
    )

    return filts, filtering_stats


def _prepare_short_reads(
    filtered_dir: Path,
    filtered_dir_rev: Path | None,
    unfilts: StrVector,
    unfilts_rev: StrVector | None,
    trunc_len: int,
    trunc_len_rev: int | None,
    trim_left: int,
    trim_left_rev: int | None,
    max_ee: float,
    max_ee_rev: float | None,
    trunc_quality: int,
    max_len: int | str,
    multithread: bool | int
) -> tuple[StrVector, StrVector | None, pd.DataFrame]:
    '''
    Filter and trim single, paired, or pyrosequencing reads.

    Parameters
    ----------
    filtered_dir : Path
        Directory in which to write filtered forward or single-end FASTQ files.
    filtered_dir_rev : Path or None
        Directory in which to write filtered reverse FASTQ files for paired
        reads.
    unfilts : StrVector
        Paths to unfiltered forward or single-end FASTQ files.
    unfilts_rev : StrVector or None
        Paths to unfiltered reverse FASTQ files for paired reads.
    trunc_len : int
        Length at which to truncate forward or single-end reads.
    trunc_len_rev : int or None
        Length at which to truncate reverse reads for paired reads.
    trim_left : int
        Number of bases to remove from the start of each forward or single-end
        read.
    trim_left_rev : int or None
        Number of bases to remove from the start of each reverse read.
    max_ee : float
        Maximum expected errors allowed in a forward or single-end read.
    max_ee_rev : float or None
        Maximum expected errors allowed in a reverse read.
    trunc_quality : int
        Quality score at which reads are truncated.
    max_len : int or str
        Maximum single-end read length, or ``"Inf"`` for no maximum.
    multithread : bool or int
        Whether to use multiple threads, or the number of threads to use.

    Returns
    -------
    filts : StrVector
        Paths to filtered forward or single-end .fastq.gz files.
    filts_rev : StrVector or None
        Paths to filtered reverse .fastq.gz files if provided.
    filtering_stats : pd.DataFrame
        Per-sample input and filtered read counts.
    '''
    filts = StrVector([
        str(filtered_dir / Path(f).name) for f in unfilts
    ])

    if unfilts_rev is not None:
        if filtered_dir_rev is None:
            raise ValueError(
                'A reverse filtered directory is required for paired reads.'
            )

        filts_rev = StrVector([
            str(filtered_dir_rev / Path(f).name) for f in unfilts_rev
        ])
        filtering_stats_r = dada2.filterAndTrim(
            unfilts, filts, unfilts_rev, filts_rev,
            truncLen=[trunc_len, trunc_len_rev],
            trimLeft=[trim_left, trim_left_rev],
            maxEE=[max_ee, max_ee_rev],
            truncQ=trunc_quality,
            rm_phix=True,
            multithread=multithread
        )
        filts_rev = StrVector([
            str(path)
            for path in _find_fastq_files(filtered_dir_rev)
        ])
    else:
        filts_rev = None
        filtering_stats_r = dada2.filterAndTrim(
            unfilts,
            filts,
            truncLen=trunc_len,
            trimLeft=trim_left,
            maxEE=max_ee,
            truncQ=trunc_quality,
            rm_phix=True,
            multithread=multithread,
            maxLen=max_len
        )

    filts = StrVector([
        str(path) for path in _find_fastq_files(filtered_dir)
    ])
    if len(filts) == 0:
        raise ValueError(
            'No reads survived filtering and trimming. '
            '(was truncLen longer than the read length?)'
        )

    filtering_stats = _robj_to_pandas_df(filtering_stats_r)
    return filts, filts_rev, filtering_stats


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


def _dereplicate_reads(filts: StrVector) -> ListVector:
    '''
    Dereplicate filtered reads into a named R list.

    Parameters
    ----------
    filts : StrVector
        Paths to filtered FASTQ files.

    Returns
    -------
    ListVector
        Named R list containing one DADA2 derep-class object per sample.
    '''
    return ListVector({
        Path(filt).name: dada2.derepFastq(filt)
        for filt in filts
    })


@dataclass(frozen=True)
class _DenoiseResults:
    '''
    Results from denoising one read direction.

    Attributes
    ----------
    samples : list[ListVector]
        Per-sample `dada-class` objects. Each `dada-class` object has the
        following key slots:
            $denoised:
                Integer vector named by inferred sequence and valued by its
                abundance.
            $map:
                Unnamed integer vector of length `derep-class$unqiues`.
                Position i is the ith unique sequence. Value at i is the index
                into `dada-class$denoised` that unique sequence i maps (was
                denoised) to. Or, NA if the unique sequence i was removed
                during denoising.

        For the remaining slots see the `dada-class` type in dada2.
    '''
    samples: list[ListVector]

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
    Dereplicate and denoise single-end, pyrosequencing, or CCS reads.

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

    if pooling_method == 'pseudo':
        dds = dada2.dada(_dereplicate_reads(filts), **kwargs)
        if 'dada' in dds.rclass:
            dds = [dds]
        else:
            dds = list(dds)
    else:
        dds = []
        for filt in filts:
            dereplicated = dada2.derepFastq(filt)
            dds.append(dada2.dada(dereplicated, **kwargs))

    return _DenoiseResults(samples=dds)


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
    filts: StrVector,
    filts_rev: StrVector,
    err: ListVector,
    err_rev: ListVector,
    pooling_method: str,
    multithread: bool | int
) -> tuple[_DenoiseResults, _DenoiseResults]:
    '''
    Dereplicate and denoise paired-end reads.

    Parameters
    ----------
    filts : StrVector
        Paths to filtered forward FASTQ files.
    filts_rev : StrVector
        Paths to filtered reverse FASTQ files.
    err : ListVector
        Learned forward-read DADA2 error model.
    err_rev : ListVector
        Learned reverse-read DADA2 error model.
    pooling_method : str
        Pooling method used during denoising. ``"pseudo"`` delegates
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
    if pooling_method == 'pseudo':
        dds_fwd = dada2.dada(
            _dereplicate_reads(filts),
            err=err,
            pool='pseudo',
            multithread=multithread,
            verbose=False
        )
        dds_rev = dada2.dada(
            _dereplicate_reads(filts_rev),
            err=err_rev,
            pool='pseudo',
            multithread=multithread,
            verbose=False
        )

        if 'dada' in dds_fwd.rclass:
            dds_fwd = [dds_fwd]
            dds_rev = [dds_rev]
        else:
            dds_fwd = list(dds_fwd)
            dds_rev = list(dds_rev)
    else:
        dds_fwd = []
        dds_rev = []
        for filt, filt_rev in zip(filts, filts_rev, strict=True):
            dds_fwd.append(dada2.dada(
                dada2.derepFastq(filt),
                err=err,
                pool=False,
                multithread=multithread,
                verbose=False
            ))
            dds_rev.append(dada2.dada(
                dada2.derepFastq(filt_rev),
                err=err_rev,
                pool=False,
                multithread=multithread,
                verbose=False
            ))

    return (
        _DenoiseResults(samples=dds_fwd),
        _DenoiseResults(samples=dds_rev)
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
    id_map = pd.DataFrame(
        {
            'temporary': pd.Series(dtype=str),
            'linked': pd.Series(dtype=str)
        }
    )

    if len(rejected) == 0:
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
    filts: StrVector,
    filts_rev: StrVector,
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
    filts : StrVector
        Paths to filtered forward FASTQ files.
    filts_rev : StrVector
        Paths to filtered reverse FASTQ files.
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

    for dd_fwd, dd_rev, filt, filt_rev in zip(
            denoised_fwd.samples, denoised_rev.samples,
            filts, filts_rev, strict=True):
        drp_fwd = dada2.derepFastq(filt)
        drp_rev = dada2.derepFastq(filt_rev)

        mp_r = dada2.mergePairs(
            dd_fwd, drp_fwd,
            dd_rev, drp_rev,
            **kwargs
        )

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
    samples: list[ListVector | RDataFrame]
) -> RMatrix:
    '''
    Construct a sequence table from per-sample DADA2 results.

    Parameters
    ----------
    samples : list[ListVector or RDataFrame]
        Per-sample DADA2 ``dada-class`` objects for single-end reads or merge
        result data frames for paired-end reads.

    Returns
    -------
    RMatrix
        Per-sample sequence table.
    '''
    return dada2.makeSequenceTable(samples)


def _validate_inputs(
    input_dir: Path,
    input_dir_rev: Path | None = None
) -> tuple[StrVector, StrVector | None]:
    '''
    Validate input directories and collect their FASTQ files.

    Parameters
    ----------
    input_dir : Path
        Directory containing forward or single-end .fastq.gz files.
    input_dir_rev : Path or None
        Optional directory containing paired reverse .fastq.gz files.

    Returns
    -------
    unfilts : StrVector
        Sorted paths to forward or single-end .fastq.gz files.
    unfilts_rev : StrVector or None
        Sorted paths to paired reverse .fastq.gz files, if provided.
    '''
    if not input_dir.exists():
        raise ValueError('Input directory does not exist.')

    unfilts = StrVector([str(path) for path in _find_fastq_files(input_dir)])

    if len(unfilts) == 0:
        raise ValueError(
            'No input files with the expected filename format (*.fastq.gz) '
            'found in forward directory.'
        )

    if input_dir_rev is not None:
        unfilts_rev = StrVector([
            str(path) for path in _find_fastq_files(input_dir_rev)
        ])

        if len(unfilts_rev) == 0:
            raise ValueError(
                'No input files with the expected filename format '
                '(*.fastq.gz) found in reverse directory.'
            )

        if len(unfilts) != len(unfilts_rev):
            raise ValueError(
                'Different numbers of forward and reverse .fastq.gz files.'
            )
    else:
        unfilts_rev = None

    return unfilts, unfilts_rev


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
    track.loc[passed_filtering, 'non-chimeric'] = list(
        non_chimeric_counts
    )

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


def _write_sequence_table(
    sequence_table: pd.DataFrame,
    filts: StrVector,
    output_path: Path
) -> None:
    '''
    Write a per-sample sequence table in the runner's TSV orientation.

    Parameters
    ----------
    sequence_table : pd.DataFrame
        Per-sample sequence table to write.
    filts : StrVector
        Filtered FASTQ paths used to derive output sample names.
    output_path : Path
        Destination for the tab-separated sequence table.

    Returns
    -------
    None
    '''
    output_table = sequence_table.T.copy()
    output_table.index.name = '#OTU ID'
    output_table.columns = [Path(filt).name for filt in filts]
    output_table.to_csv(output_path, sep='\t', index=True)


def _run_dada2(
    input_dir=None, input_dir_rev=None, output_path=None, output_track=None,
    output_err_track=None, removed_primer_dir=None, filtered_dir=None,
    filtered_dir_rev=None, forward_primer=None, reverse_primer=None,
    max_mismatch=None, indels=None, trunc_len=None, trunc_len_rev=None,
    trim_left=None, trim_left_rev=None, max_ee=None, max_ee_rev=None,
    trunc_quality=None, min_len=None, max_len=None, min_overlap=None,
    max_merge_mismatch=None, trim_overhang=None, pooling_method=None,
    chimera_method=None, min_parental_fold=None, allow_one_off=None,
    num_threads=None, learn_min_reads=None, homopolymer_gap_penalty=None,
    band_size=None, retain_unmerged=None
):
    input_dir = Path(input_dir)
    if input_dir_rev is not None:
        input_dir_rev = Path(input_dir_rev)

    unfilts, unfilts_rev = _validate_inputs(input_dir, input_dir_rev)

    if num_threads is None:
        multithread = False
    elif num_threads < 0:
        raise ValueError('Number of threads must be a positive number.')
    elif num_threads == 0:
        multithread = True
    else:
        multithread = num_threads

    if removed_primer_dir is not None:
        filts, filtering_stats = _prepare_ccs_reads(
            filtered_dir=Path(filtered_dir),
            removed_primer_dir=Path(removed_primer_dir),
            unfilts=[Path(f) for f in unfilts],
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            max_mismatch=max_mismatch,
            indels=indels,
            trunc_len=trunc_len,
            trim_left=trim_left,
            max_ee=max_ee,
            trunc_quality=trunc_quality,
            max_len=max_len,
            min_len=min_len,
            multithread=multithread
        )
        filts_rev = None
    else:
        filts, filts_rev, filtering_stats = _prepare_short_reads(
            filtered_dir=Path(filtered_dir),
            filtered_dir_rev=(
                Path(filtered_dir_rev)
                if filtered_dir_rev is not None else None
            ),
            unfilts=unfilts,
            unfilts_rev=unfilts_rev,
            trunc_len=trunc_len,
            trunc_len_rev=trunc_len_rev,
            trim_left=trim_left,
            trim_left_rev=trim_left_rev,
            max_ee=max_ee,
            max_ee_rev=max_ee_rev,
            trunc_quality=trunc_quality,
            max_len=max_len,
            multithread=multithread
        )

    error_models = _learn_error_models(
        filts=filts,
        filts_rev=filts_rev,
        learn_min_reads=learn_min_reads,
        multithread=multithread,
        pacbio=removed_primer_dir is not None,
        homopolymer_gap_penalty=homopolymer_gap_penalty,
        band_size=band_size
    )

    if input_dir_rev is None:
        denoised = _denoise_single_reads(
            filts=filts,
            err=error_models.forward,
            pooling_method=pooling_method,
            learn_min_reads=learn_min_reads,
            multithread=multithread,
            homopolymer_gap_penalty=homopolymer_gap_penalty,
            band_size=band_size
        )
        sequence_table_r = _construct_sequence_table(denoised.samples)
        denoised_counts = denoised.read_counts
        merged_counts = None
        concatenated_counts = None
    else:
        denoised_fwd, denoised_rev = _denoise_paired_reads(
            filts=filts,
            filts_rev=filts_rev,
            err=error_models.forward,
            err_rev=error_models.reverse,
            pooling_method=pooling_method,
            multithread=multithread
        )
        merged = _merge_paired_reads(
            denoised_fwd=denoised_fwd,
            denoised_rev=denoised_rev,
            filts=filts,
            filts_rev=filts_rev,
            min_overlap=min_overlap,
            max_merge_mismatch=max_merge_mismatch,
            trim_overhang=trim_overhang,
            retain_unmerged=retain_unmerged
        )
        sequence_table_r = _construct_sequence_table(merged.merged_reads)
        denoised_counts = denoised_fwd.read_counts
        merged_counts = merged.merged_counts
        concatenated_counts = (
            merged.concatenated_counts if retain_unmerged else None
        )

    sequence_table_r = _remove_chimeras(
        sequence_table=sequence_table_r,
        chimera_method=chimera_method,
        min_parental_fold=min_parental_fold,
        allow_one_off=allow_one_off,
        multithread=multithread
    )
    sequence_table = _robj_to_pandas_df(sequence_table_r)

    if input_dir_rev is not None:
        sequence_table = _restore_linked_sequences(
            sequence_table=sequence_table,
            unmerged_id_map=merged.unmerged_id_map
        )

    track = _construct_stats_table(
        filtering_stats=filtering_stats,
        denoised_counts=denoised_counts,
        non_chimeric_counts=sequence_table.sum(
            axis='columns'
        ).values,
        primer_removed=removed_primer_dir is not None,
        merged_counts=merged_counts,
        concatenated_counts=concatenated_counts
    )
    track.to_csv(output_track, sep='\t', index=True)

    error_models.stats.to_csv(output_err_track, sep='\t', index=True)
    _write_sequence_table(
        sequence_table=sequence_table,
        filts=filts,
        output_path=Path(output_path)
    )
