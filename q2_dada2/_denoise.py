# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import tempfile
import subprocess

import biom
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt)


def denoise(demultiplexed_fastq_dir: SingleLanePerSampleSingleEndFastqDirFmt,
            trunc_len: int, trim_left:int) -> biom.Table:
    with tempfile.TemporaryDirectory() as temp_dir_name:
        print(temp_dir_name)
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        cmd = ['run_dada.R', str(demultiplexed_fastq_dir), biom_fp,
               str(trunc_len), str(trim_left), temp_dir_name]
        subprocess.run(cmd, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        table = biom.Table.from_tsv(open(biom_fp), None, None, None)
    return table
