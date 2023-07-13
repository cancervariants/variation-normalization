"""A module for the Genomic Duplication Ambiguous Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, Nomenclature, GenomicDuplicationAmbiguousClassification
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier
from variation.utils import get_ambiguous_type


class GenomicDuplicationAmbiguousClassifier(Classifier):
    """The Genomic Duplication Ambiguous Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Duplication Ambiguous classification type."""
        return ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.GENOMIC_DUPLICATION_AMBIGUOUS]
        ]

    def match(self, tokens: List[Token]) -> GenomicDuplicationAmbiguousClassification:
        gene_token, genomic_dup_token = tokens
        ambiguous_type = get_ambiguous_type(
            genomic_dup_token.pos0,
            genomic_dup_token.pos1,
            genomic_dup_token.pos2,
            genomic_dup_token.pos3,
            genomic_dup_token.ambiguous_regex_type
        )

        return GenomicDuplicationAmbiguousClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=genomic_dup_token.pos0,
            pos1=genomic_dup_token.pos1,
            pos2=genomic_dup_token.pos2,
            pos3=genomic_dup_token.pos3,
            ambiguous_type=ambiguous_type
        )
