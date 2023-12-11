# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_dada2 import DADA2StatsFormat


class TestStatsBoilerplate(TestPluginBase):
    package = 'q2_dada2.tests'

    def test_dada2_stats_format_validate_positive(self):
        filenames = ['single-default-stats.tsv', 'single-override-stats.tsv',
                     'underscore-samples-stats.tsv', 'pyro-default-stats.tsv',
                     'paired-default-stats.tsv', 'paired-override-stats.tsv']
        filepaths = [self.get_data_path(os.path.join('expected', filename))
                     for filename in filenames]

        for filepath in filepaths:
            format = DADA2StatsFormat(filepath, mode='r')
            # Should pass without error
            format.validate()
            self.assertTrue(True)

    def test_dada2_stats_format_to_metadata(self):
        _, obs = self.transform_format(DADA2StatsFormat, qiime2.Metadata,
                                       os.path.join('expected',
                                                    'stats-format.tsv'))

        index = pd.Index(['L1S208', 'L1S257'], name='sample-id', dtype=object)
        cols = ['input', 'filtered', 'denoised', 'non-chimeric']
        exp_df = pd.DataFrame([[100, 99, 99, 99], [100, 98, 98, 98]],
                              index=index, columns=cols, dtype=int)
        exp = qiime2.Metadata(exp_df)
        self.assertEqual(exp, obs)

    def test_metadata_to_dada2_stats_format(self):
        transformer = self.get_transformer(qiime2.Metadata, DADA2StatsFormat)
        index = pd.Index(['L1S208', 'L1S257'], name='sample-id', dtype=object)
        cols = ['input', 'filtered', 'denoised', 'non-chimeric']
        md = qiime2.Metadata(pd.DataFrame([[100, 99, 99, 99],
                                           [100, 98, 98, 98]],
                                          index=index, columns=cols,
                                          dtype=int))
        # It shouldn't error
        transformer(md)
        self.assertTrue(True)
