# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import skbio
import biom
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)

from q2_dada2 import denoise_single, denoise_paired


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

        table, rep_seqs = denoise_single(self.demux_seqs, 100)

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))

    def test_override(self):
        with open(self.get_data_path('expected/single-override.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/single-override.fasta'),
                          'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']

        # NOTE: the test data isn't interesting enough to be impacted by
        # min_fold_parent_over_abundance.
        table, rep_seqs = denoise_single(
            self.demux_seqs, 100, trim_left=10, max_ee=10.5, trunc_q=1,
            n_threads=0, n_reads_learn=2, hashed_feature_ids=False,
            chimera_method='consensus', min_fold_parent_over_abundance=1.1)

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))

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

        # Historical NOTE: default used to be `pooled`, so the data still
        # expects that. Since this is only testing underscores, it shouldn't
        # matter much and serves as a regression test to boot.
        table, rep_seqs = denoise_single(self.demux_seqs, 100,
                                         chimera_method='pooled')

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))


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

        # NOTE: changing the chimera_method parameter doesn't impact the
        # results for this dataset
        table, rep_seqs = denoise_paired(self.demux_seqs, 150, 150)

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))

    def test_override(self):
        with open(self.get_data_path('expected/paired-override.tsv')) as fh:
            exp_table = biom.Table.from_tsv(fh, None, None, lambda x: x)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/paired-override.fasta'),
                          'fasta', constructor=skbio.DNA))
        for seq in exp_rep_seqs:
            del seq.metadata['description']

        # NOTE: the test data isn't interesting enough to be impacted by
        # chimera_method or min_fold_parent_over_abundance.
        table, rep_seqs = denoise_paired(
            self.demux_seqs, 150, 150, trim_left_f=10, trim_left_r=10,
            max_ee=20.5, trunc_q=0, n_threads=0, n_reads_learn=2,
            hashed_feature_ids=False, chimera_method='consensus',
            min_fold_parent_over_abundance=1.1)

        self.assertEqual(table, exp_table)
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))

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


if __name__ == '__main__':
    unittest.main()
