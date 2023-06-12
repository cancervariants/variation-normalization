"""A module for the DNA Coding Reference Agree Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, CdnaReferenceAgreeClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class CodingDNAReferenceAgreeClassifier(Classifier):
    """The Coding DNA Reference Agree Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA Reference Agree classification type."""
        return ClassificationType.CODING_DNA_REFERENCE_AGREE

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.CODING_DNA_REFERENCE_AGREE]
        ]

    def match(self, tokens: List[Token]) -> CdnaReferenceAgreeClassification:
        gene_token, cdna_ref_agree_token = tokens

        return CdnaReferenceAgreeClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene=gene_token,
            pos=cdna_ref_agree_token.pos
        )
