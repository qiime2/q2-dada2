# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from setuptools import setup, find_packages
import re
import ast

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('q2_dada2/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

setup(
    name="q2-dada2",
    version=version,
    packages=find_packages(),
    install_requires=['qiime >= 2.0.5', 'q2-types >= 0.0.5',
                      'biom-format >= 2.1.5, < 2.2.0'],
    author="Greg Caporaso and Benjamin Callahan",
    author_email="gregcaporaso@gmail.com",
    description="Apply DADA2 to generate denoised sequence variants. ",
    scripts=['q2_dada2/assets/profile_quality.R',
             'q2_dada2/assets/run_dada.R'],
    entry_points={
        "qiime.plugins":
        ["q2-dada2=q2_dada2.plugin_setup:plugin"]
    }
)
