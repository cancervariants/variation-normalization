"""A module for the Genomic DelIns Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, GenomicDelInsClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class GenomicDelInsClassifier(Classifier):
    """The Genomic DelIns Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic DelIns classification type."""
        return ClassificationType.GENOMIC_DELINS

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.GENOMIC_DELINS]
        ]

    def match(self, tokens: List[Token]) -> GenomicDelInsClassification:
        gene_token, genomic_delins_token = tokens

        return GenomicDelInsClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene=gene_token,
            pos0=genomic_delins_token.pos0,
            pos1=genomic_delins_token.pos1,
            inserted_sequence=genomic_delins_token.inserted_sequence
        )
