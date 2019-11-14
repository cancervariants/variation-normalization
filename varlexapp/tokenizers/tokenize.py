import re

from .amplification import Amplification
from .deletion import Deletion
from .exon import Exon
from .expression import Expression
from .fusion import Fusion
from .gain_of_function import GainOfFunction
from .gene_pair import GenePair
from .gene_symbol import GeneSymbol
from .loss_of_function import LossOfFunction
from .overexpression import OverExpression
from .protein_alternate import ProteinAlternate
from .protein_frameshift import ProteinFrameshift
from .protein_substitution import ProteinSubstitution
from .protein_termination import ProteinTermination
from .underexpression import UnderExpression
from .wild_type import WildType
from ..models import Token


from .caches import GeneSymbolCache
from .caches import AminoAcidCache

class Tokenize:
    def __init__(self, gene_file_path):
        gene_cache = GeneSymbolCache(gene_file_path)
        amino_acid_cache = AminoAcidCache()

        self.search_term_splitter = re.compile(r'\s+')

        self.tokenizers = (
                Amplification(),
                Deletion(),
                Exon(),
                Expression(),
                Fusion(),
                GainOfFunction(),
                GenePair(gene_cache),
                GeneSymbol(gene_cache),
                LossOfFunction(),
                OverExpression(),
                ProteinAlternate(amino_acid_cache),
                ProteinFrameshift(amino_acid_cache),
                ProteinSubstitution(amino_acid_cache),
                ProteinTermination(amino_acid_cache),
                UnderExpression(),
                WildType()
        )


    def perform(self, search_string):
        terms = self.search_term_splitter.split(search_string)
        tokens = list()

        for term in terms:
            matched = False
            for tokenizer in self.tokenizers:
                res = tokenizer.match(term)
                if res:
                    tokens.append(res)
                    matched = True
                    break
                else:
                    continue
            if not matched:
                tokens.append(Token('', 'Unknown', term))
        return tokens
