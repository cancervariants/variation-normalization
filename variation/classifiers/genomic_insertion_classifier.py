"""A module for the Genomic Insertion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, GenomicInsertionClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class GenomicInsertionClassifier(Classifier):
    """The Genomic Insertion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Insertion classification type."""
        return ClassificationType.GENOMIC_INSERTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.GENOMIC_INSERTION]
        ]

    def match(self, tokens: List[Token]) -> GenomicInsertionClassification:
        gene_token, genomic_ins_token = tokens

        return GenomicInsertionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=genomic_ins_token.pos0,
            pos1=genomic_ins_token.pos1,
            inserted_sequence=genomic_ins_token.inserted_sequence
        )
