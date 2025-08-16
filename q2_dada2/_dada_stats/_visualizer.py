# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
import matplotlib.pyplot as plt
import seaborn as sns
import qiime2.util
import q2templates
import qiime2

TEMPLATES = pkg_resources.resource_filename('q2_dada2._dada_stats', 'assets')


def _plot_errors(transdf, image_paths_arr, output_dir,
                 image_prefix, error_in, error_out, nominalq):

    # generates baseline graph of obs errors
    filtered_data = transdf[(transdf['from'].isin(["A", "C", "G", "T"])) &
                            (transdf['to'].isin(["A", "C", "G", "T"]))]

    g = sns.FacetGrid(filtered_data, col="Transition", col_wrap=4)
    g.map(plt.scatter, "Qual", "Observed", color="gray")

    # options for error graph; estimated line, nominal line, and input line
    if error_out is True:
        g.map(plt.plot, "Qual", "Estimated", color="gray")
    if nominalq is True:
        g.map(plt.plot, "Qual", "Nominal", color="red")
    if error_in is True:
        g.map(plt.plot, "Qual", "Input", color="blue")

    # layout options to bring it into line with the dada2 ouput version
    g.tight_layout()
    g.set(yscale="log")
    g.set_axis_labels("", "")
    g.figure.text(0.5, 0.009, 'Consensus quality score',
                  ha='center', va='center', fontsize=15)
    g.figure.text(0.009, 0.5, 'Error frequency (log10)',
                  ha='center', va='center', rotation='vertical', fontsize=15)

    # save image for html rendering
    for ext in ['png', 'svg']:
        img_fp = os.path.join(output_dir,
                              image_prefix + 'error_graph.%s' % ext)
        if ext == 'png':
            image_paths_arr.append(
                './' + image_prefix + 'error_graph.png'
            )
        plt.savefig(img_fp)
    plt.clf()
    return image_paths_arr  # pass back path to image


def plot_base_transitions(
    output_dir: str,
    base_transition_stats: qiime2.Metadata,
    nominalq: bool = False,
    error_in: bool = False,
    error_out: bool = True
):

    # intitalize indexes for html rendering
    index = os.path.join(TEMPLATES, 'index.html')

    # initalize data structures for passing to html
    image_paths_arr = []
    paired_or_not = ''
    temp_df = base_transition_stats.to_dataframe()
    # determines if the reads are paired
    if len(temp_df.columns.tolist()) > 12:
        df_fwd_subset = temp_df.iloc[:, :10]
        df_fwd_subset.columns = \
            df_fwd_subset.columns.str.replace('^F_', '', regex=True)
        df_rev_subset = temp_df.iloc[:, -10:]
        df_rev_subset.columns = \
            df_rev_subset.columns.str.replace('^R_', '', regex=True)
        image_paths_arr = _plot_errors(df_fwd_subset, image_paths_arr,
                                       output_dir, "Forward_",
                                       error_in, error_out, nominalq)
        image_paths_arr = _plot_errors(df_rev_subset, image_paths_arr,
                                       output_dir, "Reverse_",
                                       error_in, error_out, nominalq)
    else:
        image_paths_arr = _plot_errors(temp_df, image_paths_arr,
                                       output_dir, "",
                                       error_in, error_out, nominalq)

# intializes the context variable for passing to html
    context = {
        'paired_or_not': paired_or_not,
        'graph_paths': image_paths_arr
    }

    # renders the template and initalizes the
    # js sorter for the dada2 read table
    q2templates.render(index, output_dir, context=context)
