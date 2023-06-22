"""A module for the Protein Stop Gain Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, ProteinStopGainClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class ProteinStopGainClassifier(Classifier):
    """The Protein Stop Gain Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Stop Gain classification type."""
        return ClassificationType.PROTEIN_STOP_GAIN

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.PROTEIN_STOP_GAIN]
        ]

    def match(self, tokens: List[Token]):
        gene_token, protein_stop_gain_token = tokens

        return ProteinStopGainClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene=gene_token,
            pos=protein_stop_gain_token.pos,
            ref=protein_stop_gain_token.ref,
            alt=protein_stop_gain_token.alt
        )
