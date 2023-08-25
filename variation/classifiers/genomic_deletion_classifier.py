"""A module for the Genomic Deletion Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    GenomicDeletionClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType


class GenomicDeletionClassifier(Classifier):
    """The Genomic Deletion Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the genomic deletion classification.

        :return: List of list of tokens, where order matters, that represent a genomic
        deletion classification.
        """
        return [[TokenType.GENE, TokenType.GENOMIC_DELETION]]

    def match(self, tokens: List[Token]) -> GenomicDeletionClassification:
        """Return the genomic deletion classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            genomic deletion classification
        :return: genomic deletion classification for the list of matched tokens
        """
        gene_token, genomic_deletion_token = tokens

        return GenomicDeletionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=genomic_deletion_token.pos0,
            pos1=genomic_deletion_token.pos1,
            deleted_sequence=genomic_deletion_token.deleted_sequence,
        )
