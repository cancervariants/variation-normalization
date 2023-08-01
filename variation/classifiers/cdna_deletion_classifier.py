"""A module for the Cdna Deletion Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    CdnaDeletionClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType


class CdnaDeletionClassifier(Classifier):
    """The Cdna Deletion Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the cdna deletion classification.

        :return: List of list of tokens, where order matters, that represent a cdna
        deletion classification.
        """
        return [[TokenType.GENE, TokenType.CDNA_DELETION]]

    def match(self, tokens: List[Token]) -> CdnaDeletionClassification:
        """Return the cdna deletion classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            cdna deletion classification
        :return: cdna deletion classification for the list of matched tokens
        """
        gene_token, cdna_deletion_token = tokens

        return CdnaDeletionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=cdna_deletion_token.pos0,
            pos1=cdna_deletion_token.pos1,
            deleted_sequence=cdna_deletion_token.deleted_sequence,
        )
