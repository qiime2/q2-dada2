# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import skbio
import biom
import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)

from q2_dada2 import denoise_single, denoise_paired, denoise_pyro
from q2_dada2._denoise import _check_featureless_table


def _sort_seqs(seqs):
    return sorted(list(seqs), key=lambda x: x.metadata['id'])


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

        table, rep_seqs, md = denoise_single(self.demux_seqs, 100)

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(md, exp_md)

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

        # NOTE: the test data isn't interesting enough to be impacted by
        # min_fold_parent_over_abundance.
        table, rep_seqs, md = denoise_single(
            self.demux_seqs, 100, trim_left=10, max_ee=10.5, trunc_q=1,
            n_threads=1, n_reads_learn=2, hashed_feature_ids=False,
            chimera_method='consensus', min_fold_parent_over_abundance=1.1)

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(md, exp_md)

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

        # Historical NOTE: default used to be `pooled`, so the data still
        # expects that. Since this is only testing underscores, it shouldn't
        # matter much and serves as a regression test to boot.
        table, rep_seqs, md = denoise_single(self.demux_seqs, 100,
                                             chimera_method='pooled')

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(md, exp_md)

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

        table, rep_seqs, md = denoise_single(self.demux_seqs, 100,
                                             chimera_method='none')

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))

        self.assertEqual(md, exp_md)

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

        table, rep_seqs, md = denoise_single(self.demux_seqs, 100,
                                             pooling_method='pseudo')

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))

        self.assertEqual(md, exp_md)


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

        # NOTE: changing the chimera_method parameter doesn't impact the
        # results for this dataset
        table, rep_seqs, md = denoise_paired(self.demux_seqs, 150, 150)

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(md, exp_md)

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

        # NOTE: the test data isn't interesting enough to be impacted by
        # chimera_method or min_fold_parent_over_abundance.
        table, rep_seqs, md = denoise_paired(
            self.demux_seqs, 150, 150, trim_left_f=10, trim_left_r=10,
            max_ee_f=20.5, max_ee_r=20.5, trunc_q=0, n_threads=1,
            n_reads_learn=2,
            hashed_feature_ids=False, chimera_method='consensus',
            min_fold_parent_over_abundance=1.1)

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(md, exp_md)

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

        table, rep_seqs, md = denoise_paired(self.demux_seqs, 150, 150,
                                             chimera_method='none')

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(md, exp_md)


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

        table, rep_seqs, md = denoise_pyro(self.demux_seqs, 100)

        self.assertEqual(
            table,
            exp_table.sort_order(table.ids('observation'), axis='observation'))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(md, exp_md)

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


if __name__ == '__main__':
    unittest.main()
