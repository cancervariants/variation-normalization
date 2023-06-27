"""A module for the Genomic Deletion Ambiguous Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, Nomenclature, GenomicDeletionAmbiguousClassification
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class GenomicDeletionAmbiguousClassifier(Classifier):
    """The Genomic Deletion Ambiguous Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Deletion Ambiguous classification type."""
        return ClassificationType.GENOMIC_DELETION_AMBIGUOUS

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.GENOMIC_DELETION_AMBIGUOUS]
        ]

    def match(self, tokens: List[Token]) -> GenomicDeletionAmbiguousClassification:
        gene_token, genomic_dup_token = tokens

        return GenomicDeletionAmbiguousClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene=gene_token,
            pos0=genomic_dup_token.pos0,
            pos1=genomic_dup_token.pos1,
            pos2=genomic_dup_token.pos2,
            pos3=genomic_dup_token.pos3
        )
