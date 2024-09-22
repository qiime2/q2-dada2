# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
import shutil
import matplotlib.pyplot as plt
import seaborn as sns
import qiime2.util
import q2templates

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


def stats_viz(output_dir: str, dada2_stats: qiime2.Metadata,
              nominalq: bool = False, error_in: bool = False,
              error_out: bool = True):

    # intitalize indexes for html rendering
    index = os.path.join(TEMPLATES, 'index.html')
    overview_template = os.path.join(TEMPLATES, 'denoise_stats.html')
    templates = [index, overview_template]
    templates.append(os.path.join(TEMPLATES, 'error_plot_stats.html'))

    # convert stats input to dictionary for parsing
    dada_stats_dict = dict(dada2_stats)

    # initalize data structures for passing to html
    denoised_reads_data = {}
    image_paths_arr = []
    paired_or_not = ''

    # iterate over keys that link to dada2 stats files
    for keyer in dada_stats_dict.keys():

        # section for dada2 read stats
        if keyer == 'Denoised_Read_Stats':
            temp_df = dada_stats_dict[keyer].to_dataframe()
            for indexer in temp_df.index:
                row_arr = temp_df.loc[indexer]
                denoised_reads_data[indexer] \
                    = {'input': row_arr['input'],
                        'filtered': row_arr['filtered'],
                        'per_of_filt': row_arr[
                            'percentage of input passed filter'],
                        'denoised': row_arr['denoised'],
                        'non_chim': row_arr['non-chimeric'],
                        'per_of_non_chim': row_arr[
                            'percentage of input non-chimeric']}
                if 'merged' in temp_df.columns:
                    paired_or_not = 'paired'
                    denoised_reads_data[indexer]['merged'] = row_arr['merged']
                    denoised_reads_data[indexer]['per_input_merged'] = \
                        row_arr['percentage of input merged']
                elif 'primer-removed' in temp_df.columns:
                    paired_or_not = 'ccs'
                    denoised_reads_data[indexer]['primer_removed'] = \
                        row_arr['primer-removed']
                    denoised_reads_data[indexer]['per_prim_removed'] = \
                        row_arr['percentage of input primer-removed']

        # section for dada2 error plotting
        elif keyer == 'Error_Plot_Stats':
            temp_df = dada_stats_dict[keyer].to_dataframe()

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
        else:
            raise ValueError("Unknown file in metadata folder: (%r)" % keyer)
    # intializes the context variable for passing to html
    context = {
        'tabs': [{'title': 'Denoised-Stats', 'url': 'denoise_stats.html'},
                 {'title': 'Error-Plot', 'url': 'error_plot_stats.html'}],
        'denoised_data': denoised_reads_data,
        'denoised_keys': denoised_reads_data.keys(),
        'paired_or_not': paired_or_not,
        'graph_paths': image_paths_arr
    }

    # renders the template and initalizes the
    # js sorter for the dada2 read table
    q2templates.render(templates, output_dir, context=context)
    js = os.path.join(
        TEMPLATES, 'js', 'tsorter.min.js')
    os.mkdir(os.path.join(output_dir, 'js'))
    shutil.copy(js, os.path.join(output_dir, 'js', 'tsorter.min.js'))
