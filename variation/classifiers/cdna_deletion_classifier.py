"""A module for the Cdna Deletion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, CdnaDeletionClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class CdnaDeletionClassifier(Classifier):
    """The Cdna Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Cdna Deletion classification type."""
        return ClassificationType.CDNA_DELETION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.CDNA_DELETION]
        ]

    def match(self, tokens: List[Token]) -> CdnaDeletionClassification:
        gene_token, cdna_deletion_token = tokens

        return CdnaDeletionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=cdna_deletion_token.pos0,
            pos1=cdna_deletion_token.pos1,
            deleted_sequence=cdna_deletion_token.deleted_sequence
        )
