"""A module for the cDNA Reference Agree Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, CdnaReferenceAgreeClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class CdnaReferenceAgreeClassifier(Classifier):
    """The Cdna Reference Agree Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Cdna Reference Agree classification type."""
        return ClassificationType.CDNA_REFERENCE_AGREE

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.CDNA_REFERENCE_AGREE],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.CDNA_REFERENCE_AGREE]  # noqa: E501
        ]

    def match(self, tokens: List[Token]) -> CdnaReferenceAgreeClassification:
        if len(tokens) == 2:
            gene_token, cdna_ref_agree_token = tokens
        else:
            gene_token, _, cdna_ref_agree_token = tokens

        return CdnaReferenceAgreeClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos=cdna_ref_agree_token.pos
        )
