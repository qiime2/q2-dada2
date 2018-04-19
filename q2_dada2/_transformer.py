import qiime2

from q2_dada2 import DADA2StatsFormat
from q2_dada2.plugin_setup import plugin


@plugin.register_transformer
def _1(ff: DADA2StatsFormat) -> qiime2.Metadata:
    return qiime2.Metadata.load(str(ff))


@plugin.register_transformer
def _2(obj: qiime2.Metadata) -> DADA2StatsFormat:
    ff = DADA2StatsFormat()
    obj.save(str(ff))
    return ff
