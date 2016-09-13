# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import subprocess

import biom
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt)


def denoise(demultiplexed_fastq_dir: SingleLanePerSampleSingleEndFastqDirFmt)\
        -> biom.Table:
    with tempfile.NamedTemporaryFile() as f:
        cmd = ['run_dada.R', str(demultiplexed_fastq_dir), f.name]
        subprocess.run(cmd, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        table = biom.Table.from_tsv(open(f.name), None, None, None)
    return table
