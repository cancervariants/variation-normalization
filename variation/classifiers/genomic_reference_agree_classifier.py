"""A module for the Genomic Reference Agree Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, GenomicReferenceAgreeClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class GenomicReferenceAgreeClassifier(Classifier):
    """The Genomic Reference Agree Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Reference Agree classification type."""
        return ClassificationType.GENOMIC_REFERENCE_AGREE

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.CODING_DNA_REFERENCE_AGREE]
        ]

    def match(self, tokens: List[Token]) -> GenomicReferenceAgreeClassification:
        gene_token, genomic_ref_agree_token = tokens

        return GenomicReferenceAgreeClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene=gene_token,
            pos=genomic_ref_agree_token.pos
        )
