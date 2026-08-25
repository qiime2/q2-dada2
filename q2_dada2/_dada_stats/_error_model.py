# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from rpy2.rinterface import NULL
from rpy2.robjects import RObject
from rpy2.robjects.packages import importr

from q2_dada2._r_utils import _robj_to_pandas_df


dada2 = importr('dada2')


def _convert_error_matrix(robject: RObject) -> pd.DataFrame:
    df = _robj_to_pandas_df(robject)
    df.columns = [int(column) for column in df.columns]
    return df


def _melt_error_matrix(df: pd.DataFrame) -> pd.DataFrame:
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


def _error_model_to_dataframe(
    learned_errors: RObject,
    nti: tuple[str, ...] = ('A', 'C', 'G', 'T'),
    nji: tuple[str, ...] = ('A', 'C', 'G', 'T')
) -> pd.DataFrame:
    '''Convert the complete result of DADA2 error learning to statistics.'''
    acgt = {'A', 'C', 'G', 'T'}
    if not all(n in acgt for n in nti) or not all(n in acgt for n in nji):
        raise ValueError('nti and ntj must be nucleotide(s): A/C/G/T.')
    if len(set(nti)) != len(nti) or len(set(nji)) != len(nji):
        raise ValueError('nti and ntj must not contain duplicates.')

    detailed_errors = dada2.getErrors(
        learned_errors, detailed=True, enforce=False
    )

    trans_obj = detailed_errors.rx2('trans')
    err_out_obj = detailed_errors.rx2('err_out')
    if trans_obj is NULL or err_out_obj is NULL:
        raise ValueError(
            'Expected the complete result of dada2::learnErrors(), including '
            '$trans and $err_out.'
        )
    trans = _convert_error_matrix(trans_obj)
    err_out = _convert_error_matrix(err_out_obj)

    obj = detailed_errors.rx2('err_in')
    if obj is not NULL:
        if obj.rclass[0] == 'list':
            obj = obj[0]
        err_in = _convert_error_matrix(obj)
    else:
        err_in = None

    if len(trans.columns) <= 1:
        raise ValueError(
            'plotErrors only supported when using quality scores in the '
            'error model (i.e. USE_QUALS=TRUE).'
        )
    trans_df = _melt_error_matrix(trans)
    trans_df.columns = ['Transition', 'Qual', 'count']

    trans_df['from'] = trans_df['Transition'].str[0]
    trans_df['to'] = trans_df['Transition'].str[2]
    trans_df['Qual'] = pd.to_numeric(trans_df['Qual'])

    total_count = trans_df.groupby(['from', 'Qual'])['count'].sum()
    trans_df['tot'] = [
        total_count.loc[(from_base, quality)]
        for from_base, quality in zip(
            trans_df['from'], trans_df['Qual']
        )
    ]
    trans_df['Observed'] = trans_df['count'] / trans_df['tot']

    trans_df['Estimated'] = [
        err_out.loc[transition, quality]
        for transition, quality in zip(
            trans_df['Transition'], trans_df['Qual']
        )
    ]

    if err_in is not None:
        trans_df['Input'] = [
            err_in.loc[transition, quality]
            for transition, quality in zip(
                trans_df['Transition'], trans_df['Qual']
            )
        ]

    transitions = ['A2A', 'C2C', 'G2G', 'T2T']
    matching_bases = trans_df['Transition'].isin(transitions)

    trans_df['Nominal'] = (1 / 3) * (10 ** -(trans_df['Qual'] / 10))
    trans_df.loc[matching_bases, 'Nominal'] = (
        1 - (10 ** -(trans_df.loc[matching_bases, 'Qual'] / 10))
    )
    trans_df.index = pd.RangeIndex(start=1, stop=len(trans_df) + 1)

    return trans_df
