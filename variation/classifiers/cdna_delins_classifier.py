"""A module for the Cdna DelIns Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, CdnaDelInsClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class CdnaDelInsClassifier(Classifier):
    """The Cdna DelIns Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Cdna DelIns classification type."""
        return ClassificationType.CDNA_DELINS

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.CDNA_DELINS]
        ]

    def match(self, tokens: List[Token]) -> CdnaDelInsClassification:
        gene_token, cdna_delins_token = tokens

        return CdnaDelInsClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=cdna_delins_token.pos0,
            pos1=cdna_delins_token.pos1,
            inserted_sequence=cdna_delins_token.inserted_sequence
        )
