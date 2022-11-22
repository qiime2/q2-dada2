# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2


demux_single_url = \
    (f'https://docs.qiime2.org/{qiime2.__release__}/data/tutorials/'
     'moving-pictures/demux.qza')
print(demux_single_url)

demux_paired_url = \
    (f'https://docs.qiime2.org/{qiime2.__release__}/data/tutorials/'
     'atacama-soils/demux-full.qza')


def denoise_single(use):
    demux_single = use.init_artifact_from_url('demux_single', demux_single_url)

    rep_seqs, table_dada2, denoise_stats = use.action(
        use.UsageAction('dada2', 'denoise_single'),
        use.UsageInputs(
            demultiplexed_seqs=demux_single,
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


def denoise_paired(use):
    demux_paired = use.init_artifact_from_url('demux_paired', demux_paired_url)

    rep_seqs, table_dada2, denoise_stats = use.action(
        use.UsageAction('dada2', 'denoise_paired'),
        use.UsageInputs(
            demultiplexed_seqs=demux_paired,
            trunc_len_f=120,
            trunc_len_r=120,
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
