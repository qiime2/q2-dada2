# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

from ._denoise import denoise_single, denoise_paired
from ._plot import plot_qualities


__version__ = pkg_resources.get_distribution('q2-dada2').version

__all__ = ['denoise_single', 'denoise_paired', 'plot_qualities']
