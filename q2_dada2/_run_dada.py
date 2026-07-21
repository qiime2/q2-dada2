from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector, IntVector, ListVector
from rpy2.rinterface import NULL
from rpy2.robjects import pandas2ri, default_converter
from rpy2.robjects.conversion import localconverter
from pathlib import Path
import pandas as pd
import numpy as np
import os

dada2 = importr('dada2')
utils = importr('utils')
base = importr('base')


def get_n(robj):
    return sum(dada2.getUniques(robj))


def _melter(df):
    df = pd.DataFrame(df)

    melted_df = (
        df.reset_index().melt(
            id_vars='index',
            var_name='Var2',
            value_name='value'
        ).rename(columns={'index': 'Var1'})
    )

    melted_df['Var1'] = melted_df['Var1'].astype(str)
    melted_df['value'] = melted_df['value'].astype(float)
    melted_df['Var2'] = pd.to_numeric(melted_df['Var2'])

    return melted_df


def _convert_robj_to_pandas(robject) -> pd.DataFrame:
    if isinstance(robject, pd.DataFrame):
        return robject

    n_rows = int(robject.dim[0])
    n_cols = int(robject.dim[1])

    if robject.rownames is NULL:
        rownames = range(n_rows)
    else:
        rownames = list(robject.rownames)

    if robject.colnames is NULL:
        colnames = range(n_cols)
    else:
        colnames = list(robject.colnames)

    arr = np.asarray(list(robject)).reshape((n_rows, n_cols), order='F')
    return pd.DataFrame(arr, index=rownames, columns=colnames)


def _convert_error_robj(robject) -> pd.DataFrame:
    df = _convert_robj_to_pandas(robject)
    df.columns = [int(x) for x in df.columns]
    return df


def _internal_plot_errors(
    dq, nti=['A', 'C', 'G', 'T'], nji=['A', 'C', 'G', 'T']
):
    dada2 = importr('dada2')

    acgt = ['A', 'C', 'G', 'T']
    if not all(n in acgt for n in nti) or not all(n in acgt for n in nji):
        raise ValueError('nti and ntj must be nucleotide(s): A/C/G/T.')
    if len(set(nti)) != len(nti) or len(set(nji)) != len(nji):
        raise ValueError('nti and ntj must not contain duplicates.')

    dq = dada2.getErrors(dq, detailed=True, enforce=False)

    obj = dq.rx2("trans")
    if obj is not NULL:
        trans = _convert_error_robj(obj)
    else:
        trans = None

    obj = dq.rx2("err_out")
    if obj is not NULL:
        err_out = _convert_error_robj(obj)
    else:
        err_out = None

    obj = dq.rx2("err_in")
    if obj is not NULL:
        if obj.rclass[0] == "list":
            obj = obj[0]
        err_in = _convert_error_robj(obj)
    else:
        err_in = None

    if trans is not None:
        if len(trans.columns) <= 1:
            raise ValueError(
                'plotErrors only supported when using quality scores in the '
                'error model (i.e. USE_QUALS=TRUE).')
        trans_df = _melter(trans)
        trans_df.columns = ["Transition", "Qual", "count"]
    elif err_out is not None:
        if len(err_out.columns) <= 1:
            raise ValueError(
                'plotErrors only supported when using quality scores in the '
                'error model (i.e. USE_QUALS=TRUE).')
        trans_df = _melter(err_out)
        trans_df.columns = ["Transition", "Qual", "count"]
    else:
        raise ValueError(
            'Non-null observed and/or estimated error rates '
            '(dq$trans or dq$err_out) must be provided.')

    trans_df['from'] = trans_df['Transition'].str[0]
    trans_df['to'] = trans_df['Transition'].str[2]
    trans_df['Qual'] = pd.to_numeric(trans_df['Qual'])

    if trans is not None:
        total_count = trans_df.groupby(['from', 'Qual'])['count'].sum()
        trans_df['tot'] = [
            total_count.loc[(f, q)]
            for f, q in zip(trans_df['from'], trans_df['Qual'])
        ]
        trans_df['Observed'] = trans_df['count'] / trans_df['tot']
    else:
        trans_df['Observed'] = None

    if err_out is not None:
        trans_df['Estimated'] = [
            err_out.loc[t, q]
            for t, q in zip(trans_df['Transition'], trans_df['Qual'])
        ]
    else:
        trans_df['Estimated'] = None

    if err_in is not None:
        trans_df['Input'] = [
            err_in.loc[t, q]
            for t, q in zip(trans_df['Transition'], trans_df['Qual'])
        ]
    else:
        err_in = None

    transitions = ["A2A", "C2C", "G2G", "T2T"]
    mask = trans_df['Transition'].isin(transitions)

    trans_df['Nominal'] = (1/3) * (10 ** -(trans_df['Qual'] / 10))
    trans_df.loc[mask, "Nominal"] = (
        1 - (10 ** -(trans_df.loc[mask, "Qual"] / 10))
    )

    return trans_df


def remove_primers(
    filtered_dir, removed_primer_dir, unfilts, forward_primer, reverse_primer,
    max_mismatch, indels, trunc_len, trim_left, max_ee, trunc_quality,
    max_len, min_len, multithread
):
    if reverse_primer is None:
        reverse_primer = NULL
    removed_primers = [os.path.join(
        removed_primer_dir, os.path.basename(f)) for f in unfilts
    ]
    if reverse_primer is not NULL:
        reverse_primer = dada2.rc(reverse_primer)

    no_primers = dada2.removePrimers(
        fn=unfilts, fout=StrVector(removed_primers),
        **{
            'primer.f': forward_primer, 'primer.r': reverse_primer,
            'max.mismatch': max_mismatch, 'allow.indels': indels,
            'orient': True, 'verbose': True
        }
    )
    removed_primers = StrVector(sorted(
        [str(p) for p in Path(removed_primer_dir).glob('*.fastq.gz')]
    ))

    if len(removed_primers) == 0:
        raise ValueError(
            'No reads passed the Removing Primers step '
            '(Did you select the right primers?)')

    filts = StrVector(sorted([
        os.path.join(filtered_dir, os.path.basename(f))
        for f in removed_primers
    ]))
    filtered_out = dada2.filterAndTrim(
        removed_primers, filts, truncLen=trunc_len, trimLeft=trim_left,
        maxEE=max_ee, truncQ=trunc_quality, rm_phix=False,
        multithread=multithread, maxLen=max_len, minLen=min_len,
        minQ=3
    )

    out = pd.concat(
        [
            _convert_robj_to_pandas(no_primers),
            _convert_robj_to_pandas(filtered_out).iloc[:, 1]
        ],
        axis=1
    )

    return [filts, out]


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
    if not os.path.exists(input_dir):
        raise ValueError('Input directory does not exist.')

    input_dir_path = Path(input_dir)

    unfilts = StrVector(
        sorted([str(p) for p in input_dir_path.glob('*.fastq.gz')])
    )

    if len(unfilts) == 0:
        raise ValueError(
            'No input files with the expected filename format (*.fastq.gz) '
            'found in forward directory.'
        )

    if input_dir_rev is not None:
        rev_dir_path = Path(input_dir_rev)
        unfilts_rev = StrVector(
            sorted([str(p) for p in rev_dir_path.glob('*.fastq.gz')])
        )

        if len(unfilts_rev) == 0:
            raise ValueError(
                'No input files with the expected filename format '
                '(*.fastq.gz) found in reverse directory.'
            )

        if len(unfilts) != len(unfilts_rev):
            raise ValueError(
                "Different numbers of forward and reverse .fastq.gz files."
            )

    output_files = []
    if os.path.isdir(output_path):
        output_files.extend(
            os.path.join(output_path, f)
            for f in os.listdir(output_path)
        )
    if os.path.isdir(output_track):
        output_files.extend(
            os.path.join(output_track, f)
            for f in os.listdir(output_track)
        )
    if os.path.isdir(output_err_track):
        output_files.extend(
            os.path.join(output_err_track, f)
            for f in os.listdir(output_err_track)
        )

    for file in output_files:
        if os.path.isdir(file):
            raise ValueError(f'Output filename {file} is a directory.')
        elif os.path.exists(file):
            os.remove(file)

    linked_concat_delim = 'NNNNNNNNNN'

    if num_threads is None:
        multithread = False
    elif num_threads < 0:
        raise ValueError('Number of threads must be a positive number.')
    elif num_threads == 0:
        multithread = True
    else:
        multithread = num_threads

    if removed_primer_dir:
        filts, out = remove_primers(
            filtered_dir, removed_primer_dir, unfilts, forward_primer,
            reverse_primer, max_mismatch, indels, trunc_len, trim_left,
            max_ee, trunc_quality, max_len, min_len, multithread
        )
    else:
        filts = StrVector([
            os.path.join(filtered_dir, os.path.basename(f))
            for f in unfilts
        ])
        if input_dir_rev is not None:
            filts_rev = StrVector([
                os.path.join(filtered_dir_rev, os.path.basename(f))
                for f in unfilts_rev
            ])

            out = dada2.filterAndTrim(
                unfilts, filts, unfilts_rev, filts_rev,
                truncLen=[trunc_len, trunc_len_rev],
                trimLeft=[trim_left, trim_left_rev],
                maxEE=[max_ee, max_ee_rev], truncQ=trunc_quality,
                rm_phix=True, multithread=multithread
            )

            filts_rev = StrVector(sorted(
                [str(p) for p in Path(filtered_dir_rev).glob('*.fastq.gz')])
            )
        else:
            out = dada2.filterAndTrim(
                unfilts, filts, truncLen=trunc_len, trimLeft=trim_left,
                maxEE=max_ee, truncQ=trunc_quality, rm_phix=True,
                multithread=multithread, maxLen=max_len
            )

    out = _convert_robj_to_pandas(out)
    filts = StrVector(sorted(
        [str(p) for p in Path(filtered_dir).glob('*fastq.gz')])
    )

    if len(filts) == 0:
        raise ValueError(
            'No reads survived filtering and trimming. '
            '(was truncLen longer than the read length?)'
        )

    if removed_primer_dir is not None:
        err = dada2.learnErrors(
            filts, nreads=learn_min_reads,
            errorEstimationFunction=dada2.PacBioErrfun,
            multithread=multithread, BAND_SIZE=band_size
        )

        com_err_df = _internal_plot_errors(err)
    elif input_dir_rev is not None:
        err = dada2.learnErrors(
            filts, nreads=learn_min_reads, multithread=multithread
        )
        err_rev = dada2.learnErrors(
            filts_rev, nreads=learn_min_reads, multithread=multithread
        )

        err_plot_fwd = _internal_plot_errors(err)
        err_plot_fwd = err_plot_fwd.add_prefix('F_')
        err_plot_rev = _internal_plot_errors(err_rev)
        err_plot_rev = err_plot_rev.add_prefix('R_')
        com_err_df = pd.concat([err_plot_fwd, err_plot_rev], axis=1)
    else:
        kwargs = {'nreads': learn_min_reads, 'multithread': multithread}
        if homopolymer_gap_penalty is not None:
            kwargs['HOMOPOLYMER_GAP_PENALTY'] = homopolymer_gap_penalty
        if band_size is not None:
            kwargs['BAND_SIZE'] = band_size
        err = dada2.learnErrors(filts, **kwargs)
        com_err_df = _internal_plot_errors(err)

    unmerged_id_map = pd.DataFrame(
        {
            'temporary': pd.Series().astype(str),
            'linked': pd.Series().astype(str)
        }
    )

    if input_dir_rev is None:
        dds = []
        for filt in filts:
            kwargs = {
                'nreads': learn_min_reads, 'multithread': multithread,
                'err': err, 'verbose': False
            }
            if homopolymer_gap_penalty is not None:
                kwargs['HOMOPOLYMER_GAP_PENALTY'] = homopolymer_gap_penalty
            if band_size is not None:
                kwargs['BAND_SIZE'] = band_size
            derep = dada2.derepFastq(filt)
            dds.append(dada2.dada(derep, **kwargs))
        if pooling_method == 'pseudo':
            sequence_table_r = dada2.makeSequenceTable(dds)
            sequence_table = _convert_robj_to_pandas(sequence_table_r)
            prior_sequences = sequence_table.columns[
                (sequence_table > 0).sum(axis=0) >= 2
            ]
            with localconverter(default_converter + pandas2ri.converter):
                pseudo_priors = pandas2ri.py2rpy(prior_sequences)

            dds = []
            for filt in filts:
                kwargs = {
                    'nreads': learn_min_reads, 'multithread': multithread,
                    'err': err, 'verbose': False, 'priors': pseudo_priors
                }
                if homopolymer_gap_penalty is not None:
                    kwargs['HOMOPOLYMER_GAP_PENALTY'] = homopolymer_gap_penalty
                if band_size is not None:
                    kwargs['BAND_SIZE'] = band_size
                derep = dada2.derepFastq(filt)
                dds.append(dada2.dada(derep, **kwargs))

        sequence_table_r = dada2.makeSequenceTable(dds)
        sequence_table = _convert_robj_to_pandas(sequence_table_r)
    else:
        denoised_fwd = []
        merged_fwd = []
        concatenated_fwd = []
        dds_fwd = []
        dds_rev = []
        drp_fwd_list = []
        drp_rev_list = []
        mergers = []
        mergers_r = []

        for i in range(len(filts)):
            drp_fwd = dada2.derepFastq(filts[i])
            drp_fwd_list.append(drp_fwd)

            dds_fwd.append(dada2.dada(
                drp_fwd, err=err, multithread=multithread, verbose=False
            ))
            drp_rev = dada2.derepFastq(filts_rev[i])
            drp_rev_list.append(drp_rev)

            dds_rev.append(dada2.dada(
                drp_rev, err=err_rev, multithread=multithread, verbose=False
            ))

        if pooling_method == 'pseudo':
            sequence_table_fwd = dada2.makeSequenceTable(dds_fwd)
            sequence_table_fwd = _convert_robj_to_pandas(sequence_table_fwd)
            pseudo_priors_fwd = sequence_table_fwd.columns[
                (sequence_table_fwd > 0).sum(axis=0) >= 2].tolist()
            sequence_table_rev = dada2.makeSequenceTable(dds_rev)
            sequence_table_rev = _convert_robj_to_pandas(sequence_table_rev)
            pseudo_priors_rev = sequence_table_rev.columns[
                (sequence_table_rev > 0).sum(axis=0) >= 2].tolist()

            dds_fwd = []
            dds_rev = []
            drp_fwd_list = []
            drp_rev_list = []
            for i in range(len(filts)):
                drp_fwd = dada2.derepFastq(filts[i])
                drp_fwd_list.append(drp_fwd)

                dds_fwd.append(dada2.dada(
                    drp_fwd, err=err, priors=pseudo_priors_fwd,
                    multithread=multithread, verbose=False
                ))
                drp_rev = dada2.derepFastq(filts_rev[i])
                drp_rev_list.append(drp_rev)

                dds_rev.append(dada2.dada(
                    drp_rev, err=err_rev, priors=pseudo_priors_rev,
                    multithread=multithread, verbose=False
                ))

        kwargs = {}
        if min_overlap is not None:
            kwargs['minOverlap'] = min_overlap
        if max_merge_mismatch is not None:
            kwargs['maxMismatch'] = max_merge_mismatch
        if trim_overhang is not None:
            kwargs['trimOverhang'] = trim_overhang
        if retain_unmerged is not None:
            kwargs['returnRejects'] = retain_unmerged

        for i in range(len(filts)):
            drp_fwd = drp_fwd_list[i]
            drp_rev = drp_rev_list[i]

            mp_r = dada2.mergePairs(
                dds_fwd[i], drp_fwd, dds_rev[i], drp_rev, **kwargs
            )

            mp = pandas2ri.rpy2py(mp_r)
            if retain_unmerged:
                merged_fwd.append(mp.loc[mp['accept'], 'abundance'].sum())
                merge = mp.loc[mp['accept'], ['sequence', 'abundance']]
                mergers.append(merge)

                unmerged_j = mp.loc[
                    ~mp['accept'], ['forward', 'reverse', 'abundance']
                ]
                concatenated_fwd.append(unmerged_j['abundance'].sum())

                if len(unmerged_j) > 0:
                    sequence = dds_fwd[i].rx2('clustering').rx2('sequence')
                    unmerged_fwd = sequence.rx(
                        IntVector(unmerged_j['forward'].to_list())
                    )
                    sequence_rev = dds_rev[i].rx2('clustering').rx2('sequence')
                    unmerged_rev = dada2.rc(sequence_rev.rx(
                        IntVector(unmerged_j['reverse'].to_list())
                    ))

                    unmerged_temp_seqs = [
                        f'{f}{linked_concat_delim}{r}'
                        for f, r in zip(unmerged_fwd, unmerged_rev)
                    ]
                    unmerged_link_seqs = [
                        f'{f} {r}'
                        for f, r in zip(unmerged_fwd, unmerged_rev)
                    ]

                    unmerged_id_map = pd.concat(
                        [unmerged_id_map, pd.DataFrame(
                            {
                                'temporary': unmerged_temp_seqs,
                                'linked': unmerged_link_seqs
                            }
                        )], ignore_index=True
                    )
                    mergers[i] = pd.concat(
                        [mergers[i], pd.DataFrame(
                            {
                                'sequence': unmerged_temp_seqs,
                                'abundance': unmerged_j['abundance']
                            }
                        )], ignore_index=True
                    )
                    with localconverter(
                        default_converter + pandas2ri.converter
                    ):
                        mergers_r = ListVector(
                            {
                                str(i + 1): pandas2ri.py2rpy(merge)
                                for i, merge in enumerate(mergers)
                            }
                        )
            else:
                mergers_r.append(mp_r)
                merged_fwd.append(mp['abundance'].sum())

            denoised_fwd.append(get_n(dds_fwd[i]))

        sequence_table_r = dada2.makeSequenceTable(mergers_r)
        sequence_table = _convert_robj_to_pandas(sequence_table_r)

    if (
        chimera_method in ['pooled', 'consensus']
        and sequence_table.shape[1] > 0
    ):
        sequence_table_no_chimera = dada2.removeBimeraDenovo(
            sequence_table_r, method=chimera_method,
            minFoldParentOverAbundance=min_parental_fold,
            allowOneOff=allow_one_off, multithread=multithread
        )
        sequence_table_no_chimera = _convert_robj_to_pandas(
            sequence_table_no_chimera
        )
    else:
        sequence_table_no_chimera = sequence_table

    if len(unmerged_id_map) > 0:
        unmerged_id_map = unmerged_id_map.drop_duplicates()
        ambiguous_ids = unmerged_id_map.loc[
            unmerged_id_map['temporary'].duplicated(keep=False), 'temporary'
        ]
        if len(ambiguous_ids) > 0:
            raise ValueError(
                'Unable to uniquely map retained unmerged sequences from the '
                'temporary DADA2-compatible representation to linked '
                'sequences.'
            )

        unmerged_keep = sequence_table_no_chimera.columns.intersection(
            unmerged_id_map['temporary']
        ).tolist()

        if len(unmerged_keep) > 0:
            map = unmerged_id_map[
                unmerged_id_map['temporary'].isin(unmerged_keep)
            ].set_index('temporary')['linked'].to_dict()
            sequence_table_no_chimera = sequence_table_no_chimera.rename(
                columns=map
            )

    if input_dir_rev is None:
        if removed_primer_dir is not None:
            track = out.copy()
            track.columns = ["input", "primer-removed", "filtered"]
            track["denoised"] = 0
            track["non-chimeric"] = 0
        else:
            track = out.copy()
            track.columns = ["input", "filtered"]
            track["denoised"] = 0
            track["non-chimeric"] = 0
        passed_filtering = track['filtered'] > 0
        track.loc[passed_filtering, 'denoised'] = sequence_table.sum(
            axis='columns'
        ).values
        no_chimera = sequence_table_no_chimera.sum(axis='columns').values
        track.loc[passed_filtering, 'non-chimeric'] = no_chimera
        track.to_csv(output_track, sep='\t', index=True)
    else:
        if retain_unmerged:
            track = out.copy()
            track.columns = ["input", "filtered"]
            track["denoised"] = 0
            track["merged"] = 0
            track['concatenated'] = 0
            track["non-chimeric"] = 0
        else:
            track = out.copy()
            track.columns = ["input", "filtered"]
            track["denoised"] = 0
            track["merged"] = 0
            track["non-chimeric"] = 0
        passed_filtering = track['filtered'] > 0
        track.loc[passed_filtering, 'denoised'] = denoised_fwd
        print(type(merged_fwd))
        print(merged_fwd)
        track.loc[passed_filtering, 'merged'] = merged_fwd
        if retain_unmerged:
            track.loc[passed_filtering, 'concatenated'] = concatenated_fwd
        no_chimera = sequence_table_no_chimera.sum(axis='columns').values
        track.loc[passed_filtering, 'non-chimeric'] = no_chimera

        track.to_csv(output_track, sep='\t', index=True)

    com_err_df.to_csv(output_err_track, sep='\t', index=True)
    sequence_table_no_chimera = sequence_table_no_chimera.T
    col_names = [os.path.basename(file) for file in filts]
    sequence_table_no_chimera.index.name = '#OTU ID'
    sequence_table_no_chimera.columns = col_names
    sequence_table_no_chimera.to_csv(output_path, sep='\t', index=True)
