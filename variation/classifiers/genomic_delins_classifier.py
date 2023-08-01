"""A module for the Genomic DelIns Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    GenomicDelInsClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType


class GenomicDelInsClassifier(Classifier):
    """The Genomic DelIns Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the genomic delins classification.

        :return: List of list of tokens, where order matters, that represent a genomic
        delins classification.
        """
        return [[TokenType.GENE, TokenType.GENOMIC_DELINS]]

    def match(self, tokens: List[Token]) -> GenomicDelInsClassification:
        """Return the genomic delins classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            genomic delins classification
        :return: genomic delins classification for the list of matched tokens
        """
        gene_token, genomic_delins_token = tokens

        return GenomicDelInsClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=genomic_delins_token.pos0,
            pos1=genomic_delins_token.pos1,
            inserted_sequence=genomic_delins_token.inserted_sequence,
        )
