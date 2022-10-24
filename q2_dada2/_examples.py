# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2


demux_single_url = 'https://docs.qiime2.org/2022.8/data/tutorials' \
                   '/moving-pictures/demux.qza'


def artifact_from_url(url):
    def factory():
        import tempfile
        import requests

        data = requests.get(url)

        with tempfile.NamedTemporaryFile() as f:
            f.write(data.content)
            f.flush()
            result = qiime2.Artifact.load(f.name)

        return result
    return factory


def denoise_single(use):
    demux_single = use.init_artifact('demux_single',
                                     artifact_from_url(demux_single_url))

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
