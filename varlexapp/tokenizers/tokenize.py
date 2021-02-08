"""A module for tokenizing."""
import re
from typing import Iterable, List
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
from .hgvs import HGVS
from ..models import Token


from .caches import GeneSymbolCache
from .caches import AminoAcidCache


class Tokenize:
    """The tokenize class."""

    def __init__(self, gene_file_path: str) -> None:
        """Initialize the tokenize class.

        :param str gene_file_path: The path to the gene file
        """
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
                WildType(),
                HGVS()
        )

    def perform(self, search_string: str) -> Iterable[Token]:
        """Return an iterable of tokens for a given search string.

        :param str search_string: The input string to search on
        :return: An Iterable of Tokens
        """
        tokens: List[Token] = list()
        self._add_tokens(tokens, search_string)
        return tokens

    def _add_tokens(self, tokens, search_string):
        """Add tokens to a list for a given search string.

        :param list tokens: A list of tokens
        :param str search_string: The input string to search on
        """
        terms = self.search_term_splitter.split(search_string)
        for term in terms:
            if not term:
                continue
            matched = False
            for tokenizer in self.tokenizers:
                res = tokenizer.match(term)
                if res:
                    tokens.append(res)
                    if list(map(lambda t: t.token_type, tokens))[0] == 'HGVS':
                        if len(tokens) == 1:
                            self._add_tokens(tokens,
                                             search_string.split(':')[1])
                    matched = True
                    break
                else:
                    continue
            if not matched:
                tokens.append(Token('', 'Unknown', term))
