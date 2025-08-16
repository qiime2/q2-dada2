# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2

from q2_dada2 import DADA2StatsFormat, DADA2BaseTransitionStatsFormat
from q2_dada2.plugin_setup import plugin


@plugin.register_transformer
def _1(ff: DADA2StatsFormat) -> qiime2.Metadata:
    return qiime2.Metadata.load(str(ff))


@plugin.register_transformer
def _2(obj: qiime2.Metadata) -> DADA2StatsFormat:
    ff = DADA2StatsFormat()
    obj.save(str(ff))
    return ff


@plugin.register_transformer
def _3(ff: DADA2BaseTransitionStatsFormat) -> qiime2.Metadata:
    return qiime2.Metadata.load(str(ff))


@plugin.register_transformer
def _4(obj: qiime2.Metadata) -> DADA2BaseTransitionStatsFormat:
    ff = DADA2BaseTransitionStatsFormat()
    obj.save(str(ff))
    return ff
