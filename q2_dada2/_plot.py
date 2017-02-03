# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import random
import os.path
import subprocess


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


class _PlotQualView:
    """
    A very simple pass-through view which is made up of a single-end or
    paired-end directory format with a bool indicating if single or paired.
    """
    def __init__(self, directory_format, paired):
        self._directory_format = directory_format
        self.directory = str(directory_format)
        self.paired = paired


_plot_key_text = (
    "<div class='row'>\n <p class='alert alert-warning col-md-12'>\n"
    "Figure key: "
    "The distribution of quality scores at each position is shown as a "
    "grey-scale heat map, with dark colors corresponding to higher frequency. "
    "Lines show positional summary statistics: green is the mean, "
    "orange is the median, and the dashed orange lines are the 25th "
    "and 75th quantiles.\n </p>\n</div>\n")


def plot_qualities(output_dir: str, demultiplexed_seqs: _PlotQualView, n: int
                   ) -> None:
    if n < 1:
        raise ValueError("Argument to 'n' was %r, should be greater than"
                         " zero." % n)
    cmds = []
    index_f = open('%s/index.html' % output_dir, 'w')
    index_f.write('<html>\n<body>\n')
    index_f.write(_plot_key_text)
    if demultiplexed_seqs.paired:
        fps = glob.glob('%s/*.fastq.gz' % str(demultiplexed_seqs.directory))
        forward = []
        reverse = []
        for fp in fps:
            if 'R1_001.fastq' in fp:
                forward.append(fp)
            else:
                reverse.append(fp)
        fps = list(zip(sorted(forward), sorted(reverse)))
        for fwd, rev in random.sample(fps, n):
            cmds.append(['profile_quality.R', fwd, output_dir])
            cmds.append(['profile_quality.R', rev, output_dir])
            index_f.write(' <div style="display: inline-block">')
            fn = os.path.basename(fwd)
            index_f.write('  <div style="text-align: center"><b>%s</b></div>'
                          % fn.split('_', 1)[0])
            index_f.write('  <div style="float: left">')
            index_f.write('  (<a href="./%s.qprofile.pdf">PDF</a>)' % fn)
            index_f.write('  <img src="./%s.qprofile.png", width=500><br>\n'
                          % fn)
            index_f.write('  </div>')
            index_f.write('  <div style="float: left">')
            fn = os.path.basename(rev)
            index_f.write('  (<a href="./%s.qprofile.pdf">PDF</a>)' % fn)
            index_f.write('  <img src="./%s.qprofile.png", width=500><br>\n'
                          % fn)
            index_f.write('  </div>')
            index_f.write(' <span style="clear:both"></span></div>')
            index_f.write(' <hr>\n\n')
    else:
        fps = list(glob.glob('%s/*.fastq.gz'
                             % str(demultiplexed_seqs.directory)))
        for fp in random.sample(fps, n):
            cmds.append(['profile_quality.R', fp, output_dir])
            fn = os.path.basename(fp)
            index_f.write(' <div style="text-align: center"><b>%s</b></div>'
                          % fn.split('_', 1)[0])
            index_f.write(' (<a href="./%s.qprofile.pdf">PDF</a>)' % fn)
            index_f.write(' <img src="./%s.qprofile.png", width=1000><br>\n'
                          % fn)
            index_f.write(' <hr>\n\n')

    index_f.write('</body>\n</html>\n')
    index_f.close()
    run_commands(cmds)
    return None
