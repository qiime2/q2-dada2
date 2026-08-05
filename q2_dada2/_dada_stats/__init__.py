# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

__all__ = ['plot_base_transitions']


def __getattr__(name):
    if name == 'plot_base_transitions':
        from ._visualizer import plot_base_transitions
        return plot_base_transitions

    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
