import unittest

from varlexapp.tokenizers import GeneSymbol
from varlexapp.tokenizers.caches import GeneSymbolCache
from .tokenizer_base import TokenizerBase

class TestGenePairTokenizer(TokenizerBase, unittest.TestCase):

    #todo: don't hardcode this, inject with config
    def tokenizer_instance(self):
        gene_cache = GeneSymbolCache('varlexapp/data/gene_symbols.txt')
        return GeneSymbol(gene_cache)

    def token_type(self):
        return 'GeneSymbol'

    def fixture_name(self):
        return 'gene'
