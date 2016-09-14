# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import random
import os.path
import subprocess

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt)


def plot_qualities(
     output_dir: str, n: int,
     demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt)\
     -> None:
    index_f = open('%s/index.html' % output_dir, 'w')
    index_f.write('<html>\n<body>\n')
    fps = glob.glob('%s/*.fastq.gz' % str(demultiplexed_seqs))
    random.shuffle(fps)
    for fp in fps[:n]:
        cmd = ['profile_quality.R', fp, output_dir]
        fn = os.path.basename(fp)
        subprocess.run(cmd, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        index_f.write(' <b>%s</b>' % fn)
        index_f.write(' (<a href="./%s.qprofile.pdf">PDF</a>)' % fn)
        index_f.write(' <img src="./%s.qprofile.png", width=1000><br>\n' % fn)
        index_f.write(' <hr>\n\n')
    index_f.write('</body>\n</html>\n')
    index_f.close()
    return None
