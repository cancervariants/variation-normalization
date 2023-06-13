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
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.GENOMIC_DUPLICATION]
        ]

    def match(self, tokens: List[Token]) -> GenomicDuplicationClassification:
        gene_token, genomic_dup_token = tokens

        return GenomicDuplicationClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene=gene_token,
            pos0=genomic_dup_token.pos0,
            pos1=genomic_dup_token.pos1
        )
