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
from rpy2.robjects import default_converter, pandas2ri, RObject
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector, ListVector, Matrix as RMatrix
from rpy2.robjects.vectors import StrVector

from q2_dada2._dada_stats._error_model import _error_model_to_dataframe
from q2_dada2._r_utils import _robj_to_pandas_df

dada2 = importr('dada2')


def get_n(robj: RObject) -> int:
    '''
    Returns the total number of read counts in any object that contains or can
    be interpreted as a `uniques-vector`. See `dada2.getUniques` for details.
    '''
    return sum(dada2.getUniques(robj))


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
    out : pd.DataFrame
        Per-sample input, primer-removed, and filtered read counts.

    '''
    if reverse_primer is None:
        reverse_primer = NULL
    else:
        reverse_primer = dada2.rc(reverse_primer)

    removed_primers = [removed_primer_dir / f.name for f in unfilts]

    no_primers = dada2.removePrimers(
        fn=StrVector([str(f) for f in unfilts]),
        fout=StrVector([str(f) for f in removed_primers]),
        primer_fwd=forward_primer,
        primer_rev=reverse_primer,
        max_mismatch=max_mismatch,
        allow_indels=indels,
        orient=True,
        verbose=True
    )

    removed_primers = sorted(removed_primer_dir.glob('*.fastq.gz'))

    if len(removed_primers) == 0:
        raise ValueError(
            'No reads passed the Removing Primers step. Did you select the '
            'right primer(s)?'
        )

    filtered = [filtered_dir / f.name for f in removed_primers]
    filts = StrVector([str(f) for f in filtered])

    filtered_out = dada2.filterAndTrim(
        StrVector([str(f) for f in removed_primers]),
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
        str(path) for path in sorted(filtered_dir.glob('*.fastq.gz'))
    ])
    if len(filts) == 0:
        raise ValueError(
            'No reads survived filtering and trimming. '
            '(was truncLen longer than the read length?)'
        )

    out = pd.concat(
        [
            _robj_to_pandas_df(no_primers),
            _robj_to_pandas_df(filtered_out)['reads.out']
        ],
        axis=1
    )

    return filts, out


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
    out : pd.DataFrame
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
        out = dada2.filterAndTrim(
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
            for path in sorted(filtered_dir_rev.glob('*.fastq.gz'))
        ])
    else:
        filts_rev = None
        out = dada2.filterAndTrim(
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
        str(path) for path in sorted(filtered_dir.glob('*.fastq.gz'))
    ])
    if len(filts) == 0:
        raise ValueError(
            'No reads survived filtering and trimming. '
            '(was truncLen longer than the read length?)'
        )

    return filts, filts_rev, _robj_to_pandas_df(out)


@dataclass(frozen=True)
class _ErrorLearningResults:
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
        Paths to filtered reverse FASTQ files for paired-end reads.
    pacbio : bool
        Whether to learn the forward model with DADA2's PacBio error
        estimation function.
    homopolymer_gap_penalty : int or None
        Homopolymer gap penalty used when learning a single-end error model.
    band_size : int or None
        Band size used when learning a single-end or PacBio error model.

    Returns
    -------
    results : _ErrorLearningResults
        Learned forward and optional reverse models, plus their combined
        plotting statistics.
    '''
    if pacbio and filts_rev is not None:
        raise ValueError(
            'PacBio error learning does not accept reverse reads.'
        )

    kwargs = {
        'nreads': learn_min_reads,
        'multithread': multithread
    }
    if pacbio:
        kwargs['errorEstimationFunction'] = dada2.PacBioErrfun
        kwargs['BAND_SIZE'] = band_size
    else:
        if homopolymer_gap_penalty is not None:
            kwargs['HOMOPOLYMER_GAP_PENALTY'] = homopolymer_gap_penalty
        if band_size is not None:
            kwargs['BAND_SIZE'] = band_size

    forward = dada2.learnErrors(filts, **kwargs)
    forward_stats = _error_model_to_dataframe(forward)

    if filts_rev is None:
        reverse = None
        stats = forward_stats
    else:
        reverse = dada2.learnErrors(
            filts_rev,
            nreads=learn_min_reads,
            multithread=multithread
        )
        reverse_stats = _error_model_to_dataframe(reverse)
        stats = pd.concat(
            [
                forward_stats.add_prefix('F_'),
                reverse_stats.add_prefix('R_')
            ],
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
    dereplicated : ListVector
        Named R list containing one DADA2 derep-class object per sample.
    '''
    return ListVector({
        Path(filt).name: dada2.derepFastq(filt)
        for filt in filts
    })


def _denoise_single_reads(
    filts: StrVector,
    err: ListVector,
    pooling_method: str,
    learn_min_reads: int,
    multithread: bool | int,
    homopolymer_gap_penalty: int | None,
    band_size: int | None
) -> RMatrix:
    '''
    Denoise single-end, pyrosequencing, or CCS reads.

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
    sequence_table_r : RMatrix
        R sequence table constructed from the denoised samples.
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

    dereplicated = _dereplicate_reads(filts)
    dds = dada2.dada(dereplicated, **kwargs)
    return dada2.makeSequenceTable(dds)


@dataclass(frozen=True)
class _PairedDenoiseResults:
    forward: list[ListVector]
    reverse: list[ListVector]
    dereplicated_forward: list[ListVector]
    dereplicated_reverse: list[ListVector]

    @property
    def forward_read_counts(self) -> list[int]:
        return [get_n(dd) for dd in self.forward]


@dataclass(frozen=True)
class _RetainedUnmergedResults:
    mergers: pd.DataFrame
    concatenated_count: int
    id_map: pd.DataFrame


@dataclass(frozen=True)
class _PairedMergeResults:
    sequence_table: RMatrix
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
) -> _PairedDenoiseResults:
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
    results : _PairedDenoiseResults
        Denoised forward and reverse samples, their corresponding dereplicated
        reads, and the denoised forward-read count for each sample.
    '''
    pool = 'pseudo' if pooling_method == 'pseudo' else False
    dereplicated_fwd = _dereplicate_reads(filts)
    dereplicated_rev = _dereplicate_reads(filts_rev)
    dds_fwd = dada2.dada(
        dereplicated_fwd,
        err=err,
        pool=pool,
        multithread=multithread,
        verbose=False
    )
    dds_rev = dada2.dada(
        dereplicated_rev,
        err=err_rev,
        pool=pool,
        multithread=multithread,
        verbose=False
    )

    if 'dada' in dds_fwd.rclass:
        dds_fwd = [dds_fwd]
        dds_rev = [dds_rev]
    else:
        dds_fwd = list(dds_fwd)
        dds_rev = list(dds_rev)

    return _PairedDenoiseResults(
        forward=dds_fwd,
        reverse=dds_rev,
        dereplicated_forward=list(dereplicated_fwd),
        dereplicated_reverse=list(dereplicated_rev)
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
    denoised: _PairedDenoiseResults,
    min_overlap: int | None,
    max_merge_mismatch: int | None,
    trim_overhang: bool | None,
    retain_unmerged: bool | None
) -> _PairedMergeResults:
    '''
    Merge denoised paired-end reads and optionally retain rejected pairs.

    Parameters
    ----------
    denoised : _PairedDenoiseResults
        Denoised and dereplicated forward and reverse reads.
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
        Sequence table, per-sample merged and retained-unmerged counts, and
        the temporary-to-linked sequence mapping.
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

    for i in range(len(denoised.forward)):
        drp_fwd = denoised.dereplicated_forward[i]
        drp_rev = denoised.dereplicated_reverse[i]

        mp_r = dada2.mergePairs(
            denoised.forward[i], drp_fwd,
            denoised.reverse[i], drp_rev,
            **kwargs
        )

        mp = pandas2ri.rpy2py(mp_r)
        if retain_unmerged:
            merged_counts.append(
                int(mp.loc[mp['accept'], 'abundance'].sum())
            )
            retained = _retain_unmerged_pairs(
                merged_pairs=mp,
                denoised_fwd=denoised.forward[i],
                denoised_rev=denoised.reverse[i]
            )
            concatenated_counts.append(retained.concatenated_count)
            mergers.append(retained.mergers)
            unmerged_id_maps.append(retained.id_map)
        else:
            mergers_r.append(mp_r)
            merged_counts.append(int(mp['abundance'].sum()))

    if retain_unmerged:
        with localconverter(default_converter + pandas2ri.converter):
            mergers_r = ListVector(
                {
                    str(i + 1): pandas2ri.py2rpy(merger)
                    for i, merger in enumerate(mergers)
                }
            )

    if unmerged_id_maps:
        unmerged_id_map = pd.concat(unmerged_id_maps, ignore_index=True)
    else:
        unmerged_id_map = pd.DataFrame(
            {
                'temporary': pd.Series(dtype=str),
                'linked': pd.Series(dtype=str)
            }
        )

    sequence_table_r = dada2.makeSequenceTable(mergers_r)

    return _PairedMergeResults(
        sequence_table=sequence_table_r,
        merged_counts=merged_counts,
        concatenated_counts=concatenated_counts,
        unmerged_id_map=unmerged_id_map
    )


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

    unfilts = StrVector(
        sorted([str(p) for p in input_dir.glob('*.fastq.gz')])
    )

    if len(unfilts) == 0:
        raise ValueError(
            'No input files with the expected filename format (*.fastq.gz) '
            'found in forward directory.'
        )

    if input_dir_rev is not None:
        unfilts_rev = StrVector(
            sorted([str(p) for p in input_dir_rev.glob('*.fastq.gz')])
        )

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
) -> pd.DataFrame:
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
    non_chimeric_table : pd.DataFrame
        Per-sample sequence table after chimera removal, or the unchanged
        sequence table when chimera removal is skipped.
    '''
    sequence_table_df = _robj_to_pandas_df(sequence_table)
    if (
        chimera_method not in {'pooled', 'consensus'}
        or sequence_table_df.shape[1] == 0
    ):
        return sequence_table_df

    non_chimeric_table = dada2.removeBimeraDenovo(
        sequence_table,
        method=chimera_method,
        minFoldParentOverAbundance=min_parental_fold,
        allowOneOff=allow_one_off,
        multithread=multithread
    )
    return _robj_to_pandas_df(non_chimeric_table)


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
        filts, out = _prepare_ccs_reads(
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
        filts, filts_rev, out = _prepare_short_reads(
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
        sequence_table_r = _denoise_single_reads(
            filts=filts,
            err=error_models.forward,
            pooling_method=pooling_method,
            learn_min_reads=learn_min_reads,
            multithread=multithread,
            homopolymer_gap_penalty=homopolymer_gap_penalty,
            band_size=band_size
        )
        sequence_table = _robj_to_pandas_df(sequence_table_r)
    else:
        denoised = _denoise_paired_reads(
            filts=filts,
            filts_rev=filts_rev,
            err=error_models.forward,
            err_rev=error_models.reverse,
            pooling_method=pooling_method,
            multithread=multithread
        )
        merged = _merge_paired_reads(
            denoised=denoised,
            min_overlap=min_overlap,
            max_merge_mismatch=max_merge_mismatch,
            trim_overhang=trim_overhang,
            retain_unmerged=retain_unmerged
        )
        sequence_table_r = merged.sequence_table
        denoised_fwd = denoised.forward_read_counts
        sequence_table = _robj_to_pandas_df(sequence_table_r)

    sequence_table_no_chimera = _remove_chimeras(
        sequence_table=sequence_table_r,
        chimera_method=chimera_method,
        min_parental_fold=min_parental_fold,
        allow_one_off=allow_one_off,
        multithread=multithread
    )

    if input_dir_rev is not None:
        sequence_table_no_chimera = _restore_linked_sequences(
            sequence_table=sequence_table_no_chimera,
            unmerged_id_map=merged.unmerged_id_map
        )

    if input_dir_rev is None:
        denoised_counts = sequence_table.sum(axis='columns').values
        merged_counts = None
        concatenated_counts = None
    else:
        denoised_counts = denoised_fwd
        merged_counts = merged.merged_counts
        concatenated_counts = (
            merged.concatenated_counts if retain_unmerged else None
        )

    track = _construct_stats_table(
        filtering_stats=out,
        denoised_counts=denoised_counts,
        non_chimeric_counts=sequence_table_no_chimera.sum(
            axis='columns'
        ).values,
        primer_removed=removed_primer_dir is not None,
        merged_counts=merged_counts,
        concatenated_counts=concatenated_counts
    )
    track.to_csv(output_track, sep='\t', index=True)

    error_models.stats.to_csv(output_err_track, sep='\t', index=True)
    _write_sequence_table(
        sequence_table=sequence_table_no_chimera,
        filts=filts,
        output_path=Path(output_path)
    )
