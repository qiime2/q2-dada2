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

def run_commands(cmds, verbose=True):
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
    for cmd in cmds:
        if verbose:
            print("\nCommand:", end=' ')
            print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)


def plot_qualities(
     output_dir: str, n: int,
     demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt)\
     -> None:
    index_f = open('%s/index.html' % output_dir, 'w')
    index_f.write('<html>\n<body>\n')
    fps = glob.glob('%s/*.fastq.gz' % str(demultiplexed_seqs))
    random.shuffle(fps)
    cmds = []
    for fp in fps[:n]:
        cmd = ['profile_quality.R', fp, output_dir]
        cmds.append(cmd)
        fn = os.path.basename(fp)
        index_f.write(' <b>%s</b>' % fn)
        index_f.write(' (<a href="./%s.qprofile.pdf">PDF</a>)' % fn)
        index_f.write(' <img src="./%s.qprofile.png", width=1000><br>\n' % fn)
        index_f.write(' <hr>\n\n')
    index_f.write('</body>\n</html>\n')
    index_f.close()
    run_commands(cmds)
    return None
