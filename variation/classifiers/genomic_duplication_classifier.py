"""A module for the Genomic Duplication Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, Nomenclature, GenomicDuplicationClassification
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class GenomicDuplicationClassifier(Classifier):
    """The Genomic Duplication Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Duplication classification type."""
        return ClassificationType.GENOMIC_DUPLICATION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the genomic duplication classification.

        :return: List of list of tokens, where order matters, that represent a genomic
        duplication classification.
        """
        return [
            [TokenType.GENE, TokenType.GENOMIC_DUPLICATION]
        ]

    def match(self, tokens: List[Token]) -> GenomicDuplicationClassification:
        """Return the genomic duplication classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            genomic duplication classification
        :return: genomic duplication classification for the list of matched tokens
        """
        gene_token, genomic_dup_token = tokens

        return GenomicDuplicationClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=genomic_dup_token.pos0,
            pos1=genomic_dup_token.pos1
        )
