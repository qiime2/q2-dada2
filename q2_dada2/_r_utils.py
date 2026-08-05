# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from rpy2.rinterface import NULL
from rpy2.robjects import conversion, default_converter, pandas2ri, RObject
from rpy2.robjects.vectors import DataFrame as RDataFrame, Matrix as RMatrix


def _robj_to_pandas_df(robj: RObject | pd.DataFrame) -> pd.DataFrame:
    '''Convert an R dataframe or matrix into a pandas dataframe.'''
    if isinstance(robj, pd.DataFrame):
        return robj

    if not isinstance(robj, (RDataFrame, RMatrix)):
        raise ValueError(
            f'Expected `DataFrame` or `Matrix`, got `{type(robj)}` instead.'
        )

    with (default_converter + pandas2ri.converter).context():
        converted = conversion.get_conversion().rpy2py(robj)

    if isinstance(converted, pd.DataFrame):
        return converted

    n_rows, n_cols = converted.shape
    rownames = (
        range(n_rows) if robj.rownames is NULL else list(robj.rownames)
    )
    colnames = (
        range(n_cols) if robj.colnames is NULL else list(robj.colnames)
    )

    return pd.DataFrame(converted, index=rownames, columns=colnames)
