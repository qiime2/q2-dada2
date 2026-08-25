# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import skbio
import biom
from pandas.testing import assert_frame_equal

import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
from q2_types.feature_data import LinkedDNA

from q2_dada2 import denoise_single, denoise_paired, denoise_pyro, denoise_ccs
from q2_dada2._denoise import _check_featureless_table
from q2_dada2._dada_stats._visualizer import plot_base_transitions
from q2_dada2._run_dada import _ReadPaths, _prepare_paired_reads


def _sort_seqs(seqs):
    return sorted(list(seqs), key=lambda x: x.metadata['id'])


def _sort_table(table):
    return table.sort(axis="sample").sort(axis="observation")


def _assert_error_models_equal(actual, expected):
    assert_frame_equal(
        actual.to_dataframe().replace('', pd.NA),
        expected.to_dataframe().replace('', pd.NA)
    )


class TestExamples(TestPluginBase):
    package = 'q2_dada2.tests'

    def test_examples(self):
        self.execute_examples()


class TestReadPathPreparation(unittest.TestCase):

    def test_preserves_manifest_pairs_and_sample_ids(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            forward_dir = Path(temp_dir) / 'filtered-forward'
            reverse_dir = Path(temp_dir) / 'filtered-reverse'
            forward_dir.mkdir()
            reverse_dir.mkdir()

            # mix 'a' and 'z' to show filesystem sorting is no longer in effect
            unfiltered = _ReadPaths.from_manifest(pd.DataFrame(
                {
                    'forward': [
                        '/input/z-forward.fastq.gz',
                        '/input/a-forward.fastq.gz'
                    ],
                    'reverse': [
                        '/input/a-reverse.fastq.gz',
                        '/input/z-reverse.fastq.gz'
                    ]
                },
                index=['sample-b', 'sample-a']
            ))

            def filter_and_trim(fwd, filt, rev, filt_rev, **kwargs):
                self.assertEqual(list(fwd), list(unfiltered.forward))
                self.assertEqual(list(rev), list(unfiltered.reverse))
                Path(filt[0]).touch()
                Path(filt_rev[0]).touch()
                return pd.DataFrame(
                    {'reads.in': [10, 20], 'reads.out': [5, 0]},
                    index=[Path(path).name for path in fwd]
                )

            with patch(
                'q2_dada2._run_dada.dada2.filterAndTrim',
                side_effect=filter_and_trim
            ):
                filtered, filtering_stats = _prepare_paired_reads(
                    filtered_dir=forward_dir,
                    filtered_dir_rev=reverse_dir,
                    unfiltered=unfiltered,
                    trunc_len=100,
                    trunc_len_rev=100,
                    trim_left=0,
                    trim_left_rev=0,
                    max_ee=2.0,
                    max_ee_rev=2.0,
                    trunc_quality=2,
                    multithread=False
                )

            self.assertEqual(filtered.sample_ids, ('sample-b',))
            self.assertEqual(
                Path(filtered.forward[0]),
                forward_dir / 'z-forward.fastq.gz'
            )
            self.assertEqual(
                Path(filtered.reverse[0]),
                reverse_dir / 'a-reverse.fastq.gz'
            )
            self.assertEqual(
                list(filtering_stats.index), ['sample-b', 'sample-a']
            )


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

        table, rep_seqs, read_stats_md, error_model_md = denoise_single(
            self.demux_seqs, 100)
        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs), _sort_seqs(exp_rep_seqs))
        self.assertEqual(read_stats_md, exp_md)
        _assert_error_models_equal(error_model_md, exp_error_md)

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
        table, rep_seqs, read_stats_md, error_model_md = denoise_single(
            self.demux_seqs, 100, trim_left=10, max_ee=10.5, trunc_q=1,
            n_threads=1, n_reads_learn=2, hashed_feature_ids=False,
            chimera_method='consensus', min_fold_parent_over_abundance=1.1)

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(read_stats_md, exp_md)
        _assert_error_models_equal(error_model_md, exp_error_md)

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
        table, rep_seqs, read_stats_md, error_model_md = \
            denoise_single(self.demux_seqs, 100, chimera_method='pooled')

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(read_stats_md, exp_md)
        _assert_error_models_equal(error_model_md, exp_error_md)

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

        table, rep_seqs, read_stats_md, error_model_md = \
            denoise_single(self.demux_seqs, 100, chimera_method='none')

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(read_stats_md, exp_md)
        _assert_error_models_equal(error_model_md, exp_error_md)

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

        table, rep_seqs, read_stats_md, error_model_md = \
            denoise_single(self.demux_seqs, 100, pooling_method='pseudo')

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(read_stats_md, exp_md)
        _assert_error_models_equal(error_model_md, exp_error_md)


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
        table, rep_seqs, read_stats_md, error_model_md = \
            denoise_paired(self.demux_seqs, 150, 150)
        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(read_stats_md, exp_md)
        _assert_error_models_equal(error_model_md, exp_error_md)

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
        table, rep_seqs, read_stats_md, error_model_md = \
            denoise_paired(self.demux_seqs, 150, 150, retain_all_samples=False)
        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(read_stats_md, exp_md)
        _assert_error_models_equal(error_model_md, exp_error_md)

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
        table, rep_seqs, read_stats_md, error_model_md = denoise_paired(
            self.demux_seqs, 150, 150, trim_left_f=10, trim_left_r=10,
            max_ee_f=20.5, max_ee_r=20.5, trunc_q=0, n_threads=1,
            n_reads_learn=2,
            hashed_feature_ids=False, chimera_method='consensus',
            min_fold_parent_over_abundance=1.1)
        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(read_stats_md, exp_md)
        _assert_error_models_equal(error_model_md, exp_error_md)

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

        with self.assertRaisesRegex(ValueError, 'retain_unmerged'):
            denoise_paired(self.demux_seqs, 150, 150, retain_unmerged='foo')

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

        table, rep_seqs, read_stats_md, error_model_md = denoise_paired(
            self.demux_seqs, 150, 150, chimera_method='none')

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(read_stats_md, exp_md)
        _assert_error_models_equal(error_model_md, exp_error_md)


class TestDenoisePairedRetainUnmerged(TestPluginBase):
    package = 'q2_dada2.tests'

    def setUp(self):
        super().setUp()
        self.demux_seqs = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path('sample_seqs_paired'), 'r')

    def test_retain_unmerged(self):
        '''
        Smoke test that makes sure `retain_unmerged=True` yields linked
        sequences in the feature table and representative sequences.
        '''
        table, rep_seqs, _, _ = denoise_paired(
            self.demux_seqs, 150, 150,
            chimera_method='none',
            hashed_feature_ids=False,
            retain_unmerged=True
        )
        rep_seqs = list(rep_seqs)
        feature_ids = list(table.ids('observation'))

        self.assertGreater(len(table.ids('observation')), 0)
        self.assertGreater(len(table.ids('sample')), 0)

        self.assertTrue(any(' ' in seq for seq in feature_ids))
        self.assertTrue(any(' ' not in seq for seq in feature_ids))
        self.assertTrue(any(' ' in str(seq) for seq in rep_seqs))
        self.assertTrue(any(' ' not in str(seq) for seq in rep_seqs))
        self.assertTrue(all(type(seq) is LinkedDNA for seq in rep_seqs))

    def test_retain_unmerged_rescues_no_merge_run(self):
        '''
        Ensures that enabling `retain_unmerged` retains features where an
        equivalent run with `retain_unmerged` disbabled discards all features.
        '''
        shared_kwargs = dict(
            trunc_len_f=150,
            trunc_len_r=150,
            min_overlap=1000,
            chimera_method='none',
            hashed_feature_ids=False
        )

        with self.assertRaisesRegex(ValueError, 'No features remain'):
            denoise_paired(
                self.demux_seqs,
                retain_unmerged=False,
                **shared_kwargs
            )

        table, rep_seqs, read_stats_md, _ = denoise_paired(
            self.demux_seqs,
            retain_unmerged=True,
            **shared_kwargs
        )
        rep_seqs = list(rep_seqs)
        feature_ids = list(table.ids('observation'))

        self.assertGreater(len(feature_ids), 0)
        self.assertTrue(all(' ' in feature_id for feature_id in feature_ids))
        self.assertTrue(all(' ' in str(seq) for seq in rep_seqs))

        stats = read_stats_md.to_dataframe()
        self.assertEqual(list(stats.columns), [
            'input',
            'filtered',
            'percentage of input passed filter',
            'denoised',
            'merged',
            'percentage of input merged',
            'concatenated',
            'percentage of input concatenated',
            'non-chimeric',
            'percentage of input non-chimeric',
        ])
        self.assertEqual(stats['merged'].sum(), 0)
        self.assertGreater(stats['concatenated'].sum(), 0)
        self.assertEqual(
            stats['non-chimeric'].sum(),
            stats['concatenated'].sum()
        )

        exp_concat_pct = (
            stats['concatenated'] / stats['input'] * 100
        ).fillna(0).round(2)
        pd.testing.assert_series_equal(
            stats['percentage of input concatenated'],
            exp_concat_pct,
            check_names=False
        )

    def test_retain_unmerged_hashed_feature_ids_are_stable(self):
        '''
        Ensures that when `retain_unmerged` is enabled the feature hashes are
        deterministic and look like MD5 hexdigests.
        '''
        first_table, first_rep_seqs, _, _ = denoise_paired(
            self.demux_seqs, 150, 150,
            chimera_method='none',
            hashed_feature_ids=True,
            retain_unmerged=True
        )
        second_table, second_rep_seqs, _, _ = denoise_paired(
            self.demux_seqs, 150, 150,
            chimera_method='none',
            hashed_feature_ids=True,
            retain_unmerged=True
        )

        first_ids = set(first_table.ids('observation'))
        second_ids = set(second_table.ids('observation'))
        self.assertGreater(len(first_ids), 0)
        self.assertEqual(first_ids, second_ids)
        self.assertTrue(all(' ' not in feature_id for feature_id in first_ids))
        self.assertTrue(
            all(re.fullmatch(r'[0-9a-f]{32}', feature_id)
                for feature_id in first_ids)
        )

        first_rep_ids = {seq.metadata['id'] for seq in first_rep_seqs}
        second_rep_ids = {seq.metadata['id'] for seq in second_rep_seqs}
        self.assertEqual(first_ids, first_rep_ids)
        self.assertEqual(first_rep_ids, second_rep_ids)

    def test_retain_unmerged_uses_space_delimiter(self):
        '''
        Ensures that unmerged feature sequences are represented with a single
        space delimiter and that there is no leakage of DADA2's temporary
        N-separator representation.
        '''
        table, rep_seqs, _, _ = denoise_paired(
            self.demux_seqs, 150, 150,
            chimera_method='none',
            hashed_feature_ids=False,
            retain_unmerged=True
        )
        rep_seqs = list(rep_seqs)
        linked_feature_ids = [
            feature_id for feature_id in table.ids('observation')
            if ' ' in feature_id
        ]
        linked_rep_seq_strings = [
            str(seq) for seq in rep_seqs if ' ' in str(seq)
        ]

        self.assertGreater(len(linked_feature_ids), 0)
        self.assertTrue(all(
            feature_id.count(' ') == 1 for feature_id in linked_feature_ids
        ))
        self.assertTrue(all(
            'NNNNNNNNNN' not in feature_id for feature_id in linked_feature_ids
        ))
        self.assertEqual(set(linked_feature_ids), set(linked_rep_seq_strings))
        self.assertTrue(all(
            'NNNNNNNNNN' not in seq for seq in linked_rep_seq_strings
        ))

    def test_retain_unmerged_preserves_merged_features(self):
        '''
        Ensures that the set of merged features obtained is equivalent whether
        `retain_unmerged` is enabled or not.
        '''
        merged_only_table, merged_only_rep_seqs, _, _ = denoise_paired(
            self.demux_seqs, 150, 150,
            chimera_method='consensus',
            hashed_feature_ids=False,
            retain_unmerged=False
        )
        retained_table, retained_rep_seqs, _, _ = denoise_paired(
            self.demux_seqs, 150, 150,
            chimera_method='consensus',
            hashed_feature_ids=False,
            retain_unmerged=True
        )

        merged_only_ids = set(merged_only_table.ids('observation'))
        retained_merged_ids = {
            f for f in retained_table.ids('observation') if ' ' not in f
        }

        self.assertEqual(merged_only_ids, retained_merged_ids)

        merged_only_df = merged_only_table.to_dataframe(dense=True).loc[
            sorted(merged_only_ids)
        ]
        retained_merged_df = retained_table.to_dataframe(dense=True).loc[
            sorted(retained_merged_ids)
        ]
        pd.testing.assert_frame_equal(
            merged_only_df.sort_index(axis=1),
            retained_merged_df.sort_index(axis=1),
            check_dtype=False
        )

        merged_only_rep_seqs = {
            seq.metadata['id']: str(seq) for seq in merged_only_rep_seqs
        }
        retained_merged_rep_seqs = {
            seq.metadata['id']: str(seq) for seq in retained_rep_seqs
            if ' ' not in seq.metadata['id']
        }
        self.assertEqual(merged_only_rep_seqs, retained_merged_rep_seqs)

    def test_chimera_filtering_applied_to_retained_unmerged_seqs(self):
        '''
        Ensures that the rescued unmerged sequences are still processed
        by the chimera filtering algorithm, by checking that fewer unmerged
        sequences are retained when performing chimera filtering than when
        chimera filtering is not performed.
        '''
        chimera_none_table, _, chimera_none_stats, _ = denoise_paired(
            self.demux_seqs, 150, 150, min_overlap=1000, chimera_method='none',
            hashed_feature_ids=False, retain_unmerged=True
        )
        chimera_consensus_table, _, chimera_consensus_stats, _ = \
            denoise_paired(
                self.demux_seqs, 150, 150, min_overlap=1000,
                chimera_method='consensus', hashed_feature_ids=False,
                retain_unmerged=True
            )

        chimera_none_unmerged_features = {
            f for f in chimera_none_table.ids('observation') if ' ' in f
        }
        chimera_consensus_unmerged_features = {
            f for f in chimera_consensus_table.ids('observation') if ' ' in f
        }

        self.assertGreater(len(chimera_none_unmerged_features), 0)
        self.assertGreater(len(chimera_consensus_unmerged_features), 0)
        self.assertGreater(
            len(chimera_none_unmerged_features),
            len(chimera_consensus_unmerged_features)
        )
        self.assertTrue(
            chimera_consensus_unmerged_features.issubset(
                chimera_none_unmerged_features
            )
        )
        self.assertGreater(
            chimera_none_stats.to_dataframe()['non-chimeric'].sum(),
            chimera_consensus_stats.to_dataframe()['non-chimeric'].sum()
        )

        chimera_none_stats = chimera_none_stats.to_dataframe()
        chimera_consensus_stats = chimera_consensus_stats.to_dataframe()
        self.assertEqual(
            chimera_none_stats['concatenated'].sum(),
            chimera_consensus_stats['concatenated'].sum()
        )


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

        table, rep_seqs, read_stats_md, error_model_md = denoise_pyro(
            self.demux_seqs, 100)

        self.assertEqual(
            table,
            exp_table.sort_order(table.ids('observation'), axis='observation'))
        self.assertEqual(_sort_seqs(rep_seqs),
                         _sort_seqs(exp_rep_seqs))
        self.assertEqual(read_stats_md, exp_md)
        _assert_error_models_equal(error_model_md, exp_error_md)

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

        table, rep_seqs, read_stats_md, error_model_md = denoise_ccs(
            self.demux_seqs, front="AGRGTTYGATYMTGGCTCAG"
        )

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs), _sort_seqs(exp_rep_seqs))
        assert_frame_equal(
            read_stats_md.to_dataframe().sort_index(),
            exp_md.to_dataframe().sort_index()
        )
        _assert_error_models_equal(error_model_md, exp_error_md)

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

        table, rep_seqs, md, error_md = denoise_ccs(
            self.demux_seqs,
            front="AGRGTTYGATYMTGGCTCAG",
            adapter="RGYTACCTTGTTACGACTT"
        )

        self.assertEqual(_sort_table(table), _sort_table(exp_table))
        self.assertEqual(_sort_seqs(rep_seqs), _sort_seqs(exp_rep_seqs))
        read_stats_md = md
        error_model_md = error_md
        assert_frame_equal(
            read_stats_md.to_dataframe().sort_index(),
            exp_md.to_dataframe().sort_index()
        )
        _assert_error_models_equal(error_model_md, exp_error_md)


class TestVizualization(TestPluginBase):
    package = 'q2_dada2.tests'

    def setUp(self):
        super().setUp()
        self.stats_table = qiime2.Metadata.load(
                self.get_data_path('expected/single-default-error-stats.tsv'))

        self.paired_stats_table = qiime2.Metadata.load(
                self.get_data_path('expected/paired-default-error-stats.tsv'))

        self.output_dir_obj = tempfile.TemporaryDirectory(
            prefix='q2-dada2-stats-test-temp-')
        self.output_dir = self.output_dir_obj.name

    def tearDown(self):
        self.output_dir_obj.cleanup()

    def assertStat_Viz_Basics(self, viz_dir, single_or_paired):
        index_fp = os.path.join(viz_dir, 'index.html')
        self.assertTrue(os.path.exists(index_fp))
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
        plot_base_transitions(
            output_dir=self.output_dir, base_transition_stats=self.stats_table
        )
        self.assertStat_Viz_Basics(self.output_dir, True)

    def test_paired_defaults(self):
        plot_base_transitions(
            output_dir=self.output_dir,
            base_transition_stats=self.paired_stats_table
        )
        self.assertStat_Viz_Basics(self.output_dir, False)


if __name__ == '__main__':
    unittest.main()
