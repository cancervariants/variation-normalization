"""A module for the Reference Agree Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, ProteinReferenceAgreeClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class ProteinReferenceAgreeClassifier(Classifier):
    """The Reference Agree Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Reference Agree classification type."""
        return ClassificationType.PROTEIN_REFERENCE_AGREE

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.PROTEIN_REFERENCE_AGREE]
        ]

    def match(self, tokens: List[Token]):
        gene_token, protein_ref_agree_token = tokens

        return ProteinReferenceAgreeClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos=protein_ref_agree_token.pos,
            ref=protein_ref_agree_token.ref
        )
