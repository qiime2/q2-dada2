# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources

import qiime2


def _get_data_from_tests(path):
    return pkg_resources.resource_filename('q2_dada2.tests',
                                           os.path.join('data', path))


def demux_single_factory():
    return qiime2.Artifact.load(
        _get_data_from_tests('demux-01.qza'))


def denoise_single(use):
    demux = use.init_artifact('demux', demux_single_factory)

    rep_seqs, table_dada2, denoise_stats = use.action(
        use.UsageAction('dada2', 'denoise_single'),
        use.UsageInputs(
            demultiplexed_seqs=demux,
            trim_left=0,
            trunc_len=120
        ),
        use.UsageOutputNames(
            representative_sequences='representative_sequences',
            table='table',
            denoising_stats='denoising_stats'
        )
    )

    rep_seqs.assert_output_type('FeatureData[Sequence]')
    table_dada2.assert_output_type('FeatureTable[Frequency]')
    denoise_stats.assert_output_type('SampleData[DADA2Stats]')
