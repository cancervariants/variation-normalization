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
from .protein_termination import ProteinTermination
from .underexpression import UnderExpression
from .amino_acid_substitution import AminoAcidSubstitution
from .polypeptide_truncation import PolypeptideTruncation
from .silent_mutation import SilentMutation
from .wild_type import WildType
from .hgvs import HGVS
from .reference_sequence import ReferenceSequence
from variant.schemas.token_response_schema import Token, TokenMatchType


from .caches import GeneSymbolCache
from .caches import AminoAcidCache
from variant import GENE_SYMBOL_PATH


class Tokenize:
    """The tokenize class."""

    def __init__(self, gene_file_path=GENE_SYMBOL_PATH) -> None:
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
            AminoAcidSubstitution(amino_acid_cache),
            PolypeptideTruncation(amino_acid_cache),
            SilentMutation(amino_acid_cache),
            ProteinTermination(amino_acid_cache),
            UnderExpression(),
            WildType(),
            HGVS(),
            ReferenceSequence()
        )

    def perform(self, search_string: str) -> Iterable[Token]:
        """Return an iterable of tokens for a given search string.

        :param str search_string: The input string to search on
        :return: An Iterable of Tokens
        """
        tokens: List[Token] = list()
        terms = self.search_term_splitter.split(search_string)
        self._add_tokens(tokens, terms, search_string)

        # If reference sequence: Check description
        if list(map(lambda t: t.token_type, tokens)) == ['ReferenceSequence']:
            self._add_tokens(tokens,
                             [search_string.split(':')[1]],
                             search_string)
        return tokens

    def _add_tokens(self, tokens, terms, search_string):
        """Add tokens to a list for a given search string.

        :param list tokens: A list of tokens
        :param str search_string: The input string to search on
        """
        for term in terms:
            if not term:
                continue
            matched = False
            for tokenizer in self.tokenizers:
                res = tokenizer.match(term)
                if res:
                    tokens.append(res)
                    if list(map(lambda t: t.token_type, tokens))[0] == 'HGVS':
                        # Give specific type of HGVS (i.e. protein sub)
                        if len(tokens) == 1:
                            self._add_tokens(tokens,
                                             [search_string.split(':')[1]],
                                             search_string)
                    matched = True
                    break
                else:
                    continue
            if not matched:
                tokens.append(Token(
                    token='',
                    token_type='Unknown',
                    input_string=term,
                    match_type=TokenMatchType.UNSPECIFIED
                ))
