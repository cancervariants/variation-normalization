"""A module for tokenization."""
from typing import Iterable, List

from variation.schemas.token_response_schema import Token, TokenType
from variation.tokenizers import (
    HGVS,
    CdnaDeletion,
    CdnaDelIns,
    CdnaGenomicReferenceAgree,
    CdnaInsertion,
    CdnaSubstitution,
    FreeTextCategorical,
    GeneSymbol,
    GenomicDeletion,
    GenomicDelIns,
    GenomicDuplication,
    GenomicInsertion,
    GenomicSubstitution,
    GnomadVCF,
    ProteinDeletion,
    ProteinDelIns,
    ProteinInsertion,
    ProteinReferenceAgree,
    ProteinSubstitution,
)

r"(\((\?|d+)_(\?|\d+)\))_(\((\?|\d+)_(\?|\d+)\))dup"


class Tokenize:
    """The tokenize class."""

    def __init__(self, gene_symbol: GeneSymbol) -> None:
        """Initialize the tokenize class."""
        self.gene_symbol = gene_symbol
        self.tokenizers = [
            HGVS(),
            GnomadVCF(),
            self.gene_symbol,
            FreeTextCategorical(),
            # Substitution
            ProteinSubstitution(),
            GenomicSubstitution(),
            CdnaSubstitution(),
            # Reference Agree
            ProteinReferenceAgree(),
            CdnaGenomicReferenceAgree(),
            # Delins
            ProteinDelIns(),
            CdnaDelIns(),
            GenomicDelIns(),
            # Deletion
            ProteinDeletion(),
            CdnaDeletion(),
            GenomicDeletion(),
            # Insertion
            ProteinInsertion(),
            CdnaInsertion(),
            GenomicInsertion(),
            # Duplication
            GenomicDuplication(),
        ]

    def perform(self, search_string: str, warnings: List[str]) -> Iterable[Token]:
        """Return an iterable of tokens for a given search string.

        :param str search_string: The input string to search on
        :param List warnings: List of warnings
        :return: An Iterable of Tokens
        """
        terms = search_string.split()

        tokens: List[Token] = []
        for term in terms:
            if not term:
                continue

            matched = False
            for tokenizer in self.tokenizers:
                res = tokenizer.match(term)
                if res:
                    if isinstance(res, List):
                        for r in res:
                            tokens.append(r)
                            if not matched:
                                matched = True
                    else:
                        tokens.append(res)
                        matched = True
                        break

            if not matched:
                warnings.append(f"Unable to tokenize: {term}")
                tokens.append(
                    Token(token=term, token_type=TokenType.UNKNOWN, input_string=term)
                )

        return tokens
