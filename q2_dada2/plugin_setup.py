# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime.plugin
from q2_types.per_sample_sequences import SequencesWithQuality
from q2_types.sample_data import SampleData
from q2_types.feature_table import FeatureTable, Frequency


import q2_dada2


plugin = qiime.plugin.Plugin(
    name='dada2',
    version=q2_dada2.__version__,
    website='http://benjjneb.github.io/dada2/',
    package='q2_dada2',
    user_support_text=("To get help with DADA2, post to the DADA2 issue "
                       "tracker: https://github.com/benjjneb/dada2/issues"),
    citation_text=("DADA2: High-resolution sample inference from Illumina "
                   "amplicon data. Benjamin J Callahan, Paul J McMurdie, "
                   "Michael J Rosen, Andrew W Han, Amy Jo A Johnson, "
                   "Susan P Holmes. Nature Methods 13, 581–583 (2016) "
                   "doi:10.1038/nmeth.3869.")
)


plugin.methods.register_function(
    function=q2_dada2.denoise,
    inputs={'demultiplexed_fastq_dir': SampleData[SequencesWithQuality]},
    parameters={'trunc_len': qiime.plugin.Int,
                'trim_left': qiime.plugin.Int},
    outputs=[('table', FeatureTable[Frequency])],
    name='Denoise and dereplicate',
    description=('This method denoises sequences, dereplicates them, and '
                 'filters chimeras.')
)

plugin.visualizers.register_function(
    function=q2_dada2.plot_qualities,
    inputs={'demultiplexed_fastq_dir': SampleData[SequencesWithQuality]},
    parameters={'n': qiime.plugin.Int},
    name='Plot positional qualitites',
    description=('Plots positional quality scores for n samples selected '
                 'at random from the input data.')
)
