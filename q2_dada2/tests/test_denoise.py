# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import tempfile
import pandas as pd
import skbio
import biom
import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)

from q2_dada2 import denoise_single, denoise_paired, denoise_pyro, denoise_ccs
from q2_dada2._denoise import _check_featureless_table
from q2_dada2._dada_stats._visualizer import (stats_viz)
import os


def _sort_seqs(seqs):
    return sorted(list(seqs), key=lambda x: x.metadata['id'])


def _sort_table(table):
    return table.sort(axis="sample").sort(axis="observation")


class TestExamples(TestPluginBase):
    package = 'q2_dada2.tests'

    def test_examples(self):
        self.execute_examples()


class TestDenoiseSingle(TestPluginBase):
    package = 'q2_dada2.tests'

    def setUp(self):
        super().setUp()
        self.demux_seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('sample_seqs_single'), 'r')

    def test_defaults(self):
        with open(self.get_data_path('expected/single-default.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/single-default.fasta'),
                          'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']
        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/single-default-stats.tsv'))
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/single-default-error-stats.tsv'))

        table, rep_seqs, md = denoise_single(self.demux_seqs, 100)
        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs), _sort_seqs(exp_rep_seqs))
        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(error_model_md.to_dataframe().replace('', pd.NA, inplace=True),
                         exp_error_md.to_dataframe().replace('', pd.NA, inplace=True))

    def test_override(self):
        with open(self.get_data_path('expected/single-override.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/single-override.fasta'),
                          'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']
        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/single-override-stats.tsv'))
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/single-override-error-stats.tsv'))

        # NOTE: the test data isn't interesting enough to be impacted by
        # min_fold_parent_over_abundance.
        table, rep_seqs, md = denoise_single(
            self.demux_seqs, 100, trim_left=10, max_ee=10.5, trunc_q=1,
            n_threads=1, n_reads_learn=2, hashed_feature_ids=False,
            chimera_method='consensus', min_fold_parent_over_abundance=1.1)

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(error_model_md.to_dataframe().replace('', pd.NA, inplace=True),
                         exp_error_md.to_dataframe().replace('', pd.NA, inplace=True))

    def test_mixed_barcodes_and_ids(self):
        demux_seqs = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path('mixed_barcodes_and_ids'), 'r')

        denoise_paired(demux_seqs, 150, 150)

        self.assertTrue(True)

    def test_all_reads_filtered(self):
        with self.assertRaisesRegex(ValueError, 'filter'):
            denoise_single(self.demux_seqs, 10000)

    def test_bad_values_fail(self):
        # Just confirm that the machinery works, anything more specific is just
        # restating the _valid_inputs dict which is more declarative than a
        # unit-test anyways.
        with self.assertRaisesRegex(ValueError, 'trunc_len'):
            denoise_single(self.demux_seqs, -1)

        with self.assertRaisesRegex(ValueError, 'n_reads_learn'):
            denoise_single(self.demux_seqs, 100, n_reads_learn=0)

        with self.assertRaisesRegex(ValueError, 'consensus'):
            denoise_single(self.demux_seqs, 100, chimera_method='foo')

    def test_trim_left_bigger_than_trunc_len(self):
        with self.assertRaisesRegex(ValueError, 'trim_left'):
            denoise_single(self.demux_seqs, 100, trim_left=100)

        # Shouldn't fail when `trunc_len=0`
        denoise_single(self.demux_seqs, 0, trim_left=100)

    def test_underscore_samples(self):
        self.demux_seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('underscore_samples'), 'r')

        with open(self.get_data_path('expected/underscore-samples.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(
                self.get_data_path('expected/underscore-samples.fasta'),
                'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']
        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/underscore-samples-stats.tsv'))
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/single-default-error-stats.tsv'))

        # Historical NOTE: default used to be `pooled`, so the data still
        # expects that. Since this is only testing underscores, it shouldn't
        # matter much and serves as a regression test to boot.
        table, rep_seqs, md = denoise_single(self.demux_seqs, 100,
                                             chimera_method='pooled')

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(error_model_md.to_dataframe().replace('', pd.NA, inplace=True),
                         exp_error_md.to_dataframe().replace('', pd.NA, inplace=True))

    def test_no_chimera_method(self):
        with open(self.get_data_path('expected/single-default.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/single-default.fasta'),
                          'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']
        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/single-default-stats.tsv'))
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/single-default-error-stats.tsv'))

        table, rep_seqs, md = denoise_single(self.demux_seqs, 100,
                                             chimera_method='none')

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))

        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(error_model_md.to_dataframe().replace('', pd.NA, inplace=True),
                         exp_error_md.to_dataframe().replace('', pd.NA, inplace=True))

    def test_pseudo_pooling(self):
        with open(self.get_data_path('expected/single-pseudo.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/single-pseudo.fasta'),
                          'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']
        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/single-pseudo-stats.tsv'))
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/single-default-error-stats.tsv'))

        table, rep_seqs, md = denoise_single(self.demux_seqs, 100,
                                             pooling_method='pseudo')

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))

        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(error_model_md.to_dataframe().replace('', pd.NA, inplace=True),
                         exp_error_md.to_dataframe().replace('', pd.NA, inplace=True))


class TestDenoisePaired(TestPluginBase):
    package = 'q2_dada2.tests'

    def setUp(self):
        super().setUp()
        self.demux_seqs = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path('sample_seqs_paired'), 'r')

    def test_defaults(self):
        with open(self.get_data_path('expected/paired-default.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/paired-default.fasta'),
                          'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']
        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/paired-default-stats.tsv'))
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/paired-default-error-stats.tsv'))
        # NOTE: changing the chimera_method parameter doesn't impact the
        # results for this dataset
        table, rep_seqs, md = denoise_paired(self.demux_seqs, 150, 150)
        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(error_model_md.to_dataframe().replace('', pd.NA, inplace=True),
                         exp_error_md.to_dataframe().replace('', pd.NA, inplace=True))

    def test_remove_empty(self):
        with open(self.get_data_path('expected/paired-remove-empty-default.tsv'
                                     )) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/paired-default.fasta'),
                          'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']
        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/paired-default-stats.tsv'))
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/paired-default-error-stats.tsv'))
        # NOTE: changing the chimera_method parameter doesn't impact the
        # results for this dataset
        table, rep_seqs, md = denoise_paired(self.demux_seqs, 150, 150,
                                             retain_all_samples=False)
        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(error_model_md.to_dataframe().replace('', pd.NA, inplace=True),
                         exp_error_md.to_dataframe().replace('', pd.NA, inplace=True))

    def test_override(self):
        with open(self.get_data_path('expected/paired-override.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/paired-override.fasta'),
                          'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']
        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/paired-override-stats.tsv'))
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/paired-override-error-stats.tsv'))

        # NOTE: the test data isn't interesting enough to be impacted by
        # chimera_method or min_fold_parent_over_abundance.
        table, rep_seqs, md = denoise_paired(
            self.demux_seqs, 150, 150, trim_left_f=10, trim_left_r=10,
            max_ee_f=20.5, max_ee_r=20.5, trunc_q=0, n_threads=1,
            n_reads_learn=2,
            hashed_feature_ids=False, chimera_method='consensus',
            min_fold_parent_over_abundance=1.1)
        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(error_model_md.to_dataframe().replace('', pd.NA, inplace=True),
                         exp_error_md.to_dataframe().replace('', pd.NA, inplace=True))

    def test_all_reads_filtered(self):
        with self.assertRaisesRegex(ValueError, 'filter'):
            denoise_paired(self.demux_seqs, 10000, 10000)

        with self.assertRaisesRegex(ValueError, 'filter'):
            denoise_paired(self.demux_seqs, 150, 10000)

        with self.assertRaisesRegex(ValueError, 'filter'):
            denoise_paired(self.demux_seqs, 10000, 150)

    def test_bad_values_fail(self):
        # Just confirm that the machinery works, anything more specific is just
        # restating the _valid_inputs dict which is more declarative than a
        # unit-test anyways.
        with self.assertRaisesRegex(ValueError, 'trunc_len_f'):
            denoise_paired(self.demux_seqs, -1, 150)

        with self.assertRaisesRegex(ValueError, 'trunc_len_r'):
            denoise_paired(self.demux_seqs, 150, -1)

        with self.assertRaisesRegex(ValueError, 'n_reads_learn'):
            denoise_paired(self.demux_seqs, 150, 150, n_reads_learn=0)

        with self.assertRaisesRegex(ValueError, 'consensus'):
            denoise_single(self.demux_seqs, 150, 150, chimera_method='foo')

    def test_trim_left_bigger_than_trunc_len(self):
        with self.assertRaisesRegex(ValueError, 'trim_left_f'):
            denoise_paired(self.demux_seqs, 150, 150, trim_left_f=150)

        with self.assertRaisesRegex(ValueError, 'trim_left_r'):
            denoise_paired(self.demux_seqs, 150, 150, trim_left_r=150)

        # Shouldn't fail when `trunc_len_f=0`
        denoise_paired(self.demux_seqs, 0, 150, trim_left_f=10)
        # Shouldn't fail when `trunc_len_r=0`
        denoise_paired(self.demux_seqs, 150, 0, trim_left_r=10)

    def test_no_chimera_method(self):
        with open(self.get_data_path('expected/paired-default.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/paired-default.fasta'),
                          'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']
        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/paired-default-stats.tsv'))
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/paired-default-error-stats.tsv'))

        table, rep_seqs, md = denoise_paired(self.demux_seqs, 150, 150,
                                             chimera_method='none')

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(error_model_md.to_dataframe().replace('', pd.NA, inplace=True),
                         exp_error_md.to_dataframe().replace('', pd.NA, inplace=True))


# More thorough tests exist in TestDenoiseSingle --- denoise-pyro is basically
# just a variation of denoise-single. These tests should serve as regression
# or integration tests (depending on perspective).
class TestDenoisePyro(TestPluginBase):
    package = 'q2_dada2.tests'

    def setUp(self):
        super().setUp()
        # Reusing the single-end reads for this test suite
        self.demux_seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('sample_seqs_single'), 'r')

    def test_defaults(self):
        with open(self.get_data_path('expected/pyro-default.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/pyro-default.fasta'),
                          'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']
        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/pyro-default-stats.tsv'))
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/pyro-default-error-stats.tsv'))

        table, rep_seqs, md = denoise_pyro(self.demux_seqs, 100)

        self.assertEqual(
            table,
            exp_table.sort_order(table.ids('observation'), axis='observation'))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(error_model_md.to_dataframe().replace('', pd.NA, inplace=True),
                         exp_error_md.to_dataframe().replace('', pd.NA, inplace=True))

    def test_trunc_len_bigger_than_max_len(self):
        with self.assertRaisesRegex(ValueError, 'max_len'):
            denoise_pyro(self.demux_seqs, 100, max_len=99)

        # Shouldn't fail when max_len > trunc_len
        denoise_pyro(self.demux_seqs, 100, max_len=160)


class TestUtils(TestPluginBase):
    package = 'q2_dada2.tests'

    def test_check_featureless_table_single_feature(self):
        fp = self.get_data_path('single_feature.tsv')

        # should not raise an error
        _check_featureless_table(fp)

        self.assertTrue(True)

    def test_check_featureless_table_no_features(self):
        fp = self.get_data_path('no_asvs.tsv')

        with self.assertRaisesRegex(ValueError, "No features"):
            _check_featureless_table(fp)


class TestDenoiseCCS(TestPluginBase):
    package = 'q2_dada2.tests'

    def setUp(self):
        super().setUp()
        self.demux_seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('sample_seqs_ccs'), 'r')

    def test_default(self):
        with open(self.get_data_path('expected/ccs-default.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)

        exp_rep_seqs = list(
            skbio.io.read(
                self.get_data_path('expected/ccs-default.fasta'),
                'fasta',
                constructor=skbio.DNA
            )
        )

        for seq in exp_rep_seqs:
            del seq.metadata['description']

        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/ccs-default-stats.tsv')
        )
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/ccs-default-error-stats.tsv'))

        table, rep_seqs, md = denoise_ccs(
            self.demux_seqs, front="AGRGTTYGATYMTGGCTCAG"
        )

        self.assertEqual(
            table,
            exp_table.sort_order(
                table.ids('observation'),
                axis='observation'
            )
        )
        self.assertEqual(_sort_seqs(rep_seqs), _sort_seqs(exp_rep_seqs))
        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        df_err_md = \
            error_model_md.to_dataframe().replace('', pd.NA, inplace=True)
        df_err_exp_md = \
            exp_error_md.to_dataframe().replace('', pd.NA, inplace=True)
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(df_err_md, df_err_exp_md)

    def test_with_reverse_primer(self):
        with open(self.get_data_path('expected/ccs-reverse-primer.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)

        exp_rep_seqs = list(
            skbio.io.read(
                self.get_data_path('expected/ccs-reverse-primer.fasta'),
                'fasta',
                constructor=skbio.DNA
            )
        )

        for seq in exp_rep_seqs:
            del seq.metadata['description']

        exp_md = qiime2.Metadata.load(
            self.get_data_path('expected/ccs-reverse-primer-stats.tsv')
        )
        exp_error_md = qiime2.Metadata.load(
            self.get_data_path('expected/ccs-reverse-primer-error-stats.tsv'))

        table, rep_seqs, md = denoise_ccs(
            self.demux_seqs,
            front="AGRGTTYGATYMTGGCTCAG",
            adapter="RGYTACCTTGTTACGACTT"
        )

        self.assertEqual(
            table,
            exp_table.sort_order(
                table.ids('observation'),
                axis='observation'
            )
        )
        self.assertEqual(_sort_seqs(rep_seqs), _sort_seqs(exp_rep_seqs))
        read_stats_md = dict(md)["Denoised_Read_Stats"]
        error_model_md = dict(md)["Error_Plot_Stats"]
        df_err_md = \
            error_model_md.to_dataframe().replace('', pd.NA, inplace=True)
        df_err_exp_md = \
            exp_error_md.to_dataframe().replace('', pd.NA, inplace=True)
        self.assertEqual(read_stats_md, exp_md)
        self.assertEqual(df_err_md, df_err_exp_md)


class TestVizualization(TestPluginBase):
    package = 'q2_dada2.tests'

    def setUp(self):
        super().setUp()
        self.stats_table = {'Denoised_Read_Stats': qiime2.Metadata.load(
            self.get_data_path('expected/single-default-stats.tsv')),
            'Error_Plot_Stats': qiime2.Metadata.load(
                self.get_data_path('expected/single-default-error-stats.tsv'))}

        self.paired_stats_table = {'Denoised_Read_Stats': qiime2.Metadata.load(
            self.get_data_path('expected/paired-default-stats.tsv')),
            'Error_Plot_Stats':  qiime2.Metadata.load(
                self.get_data_path('expected/paired-default-error-stats.tsv'))}

        self.output_dir_obj = tempfile.TemporaryDirectory(
            prefix='q2-dada2-stats-test-temp-')
        self.output_dir = self.output_dir_obj.name

    def tearDown(self):
        self.output_dir_obj.cleanup()

    def assertStat_Viz_Basics(self, viz_dir, single_or_paired):
        index_fp = os.path.join(viz_dir, 'index.html')
        self.assertTrue(os.path.exists(index_fp))
        with open(index_fp, 'r') as fh:
            index_contents = fh.read()
        self.assertIn('./denoise_stats.html', index_contents)
        self.assertIn('./error_plot_stats.html', index_contents)
        if single_or_paired is True:
            self.assertTrue(
                os.path.exists(os.path.join(viz_dir, 'error_graph.png')))
        else:
            self.assertTrue(
                os.path.exists(
                    os.path.join(viz_dir, 'Reverse_error_graph.png')))
            self.assertTrue(
                os.path.exists(
                    os.path.join(viz_dir, 'Forward_error_graph.png')))

    def test_defaults(self):
        stats_viz(output_dir=self.output_dir,
                  dada2_stats=self.stats_table)
        self.assertStat_Viz_Basics(self.output_dir, True)

    def test_paired_defaults(self):
        stats_viz(output_dir=self.output_dir,
                  dada2_stats=self.paired_stats_table)
        self.assertStat_Viz_Basics(self.output_dir, False)


if __name__ == '__main__':
    unittest.main()
