from qiime2.plugin import SemanticType, model
from q2_types.sample_data import SampleData


DADA2Stats = SemanticType('DADA2Stats', variant_of=SampleData.field['type'])


class DADA2StatsFormat(model.TextFileFormat):
    def validate(*args):
        pass


DADA2StatsDirFmt = model.SingleFileDirectoryFormat(
    'DADA2StatsDirFmt', 'stats.tsv', DADA2StatsFormat)
