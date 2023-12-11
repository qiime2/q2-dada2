# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

import versioneer


setup(
    name="q2-dada2",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url="https://qiime2.org",
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Greg Caporaso and Benjamin Callahan",
    author_email="gregcaporaso@gmail.com",
    description="Apply DADA2 to generate denoised sequence variants. ",
    scripts=['q2_dada2/assets/run_dada.R'],
    package_data={
        'q2_dada2': ['citations.bib'],
        'q2_dada2.tests': ['data/*',
                           'data/expected/*',
                           'data/mixed_barcodes_and_ids/*',
                           'data/underscore_samples/*',
                           'data/sample_seqs_single/*',
                           'data/sample_seqs_ccs/*',
                           'data/sample_seqs_paired/*']
    },
    entry_points={
        "qiime2.plugins":
        ["q2-dada2=q2_dada2.plugin_setup:plugin"]
    },
    zip_safe=False,
)
