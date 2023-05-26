"""A module for the Genomic Deletion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, GenomicDeletionClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class GenomicDeletionClassifier(Classifier):
    """The Genomic Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Deletion classification type."""
        return ClassificationType.GENOMIC_DELETION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.GENOMIC_DELETION]
        ]

    def match(self, tokens: List[Token]):
        gene_token, cdna_deletion_token = tokens

        return GenomicDeletionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene=gene_token,
            pos0=cdna_deletion_token.pos0,
            pos1=cdna_deletion_token.pos1,
            deleted_sequence=cdna_deletion_token.deleted_sequence
        )
