"""A module for the Protein Stop Gain Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType,
    ProteinStopGainClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class ProteinStopGainClassifier(Classifier):
    """The Protein Stop Gain Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Stop Gain classification type."""
        return ClassificationType.PROTEIN_STOP_GAIN

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the protein stop gain classification.

        :return: List of list of tokens, where order matters, that represent a protein
        stop gain classification.
        """
        return [[TokenType.GENE, TokenType.PROTEIN_STOP_GAIN]]

    def match(self, tokens: List[Token]) -> ProteinStopGainClassification:
        """Return the protein stop gain classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            protein stop gain classification
        :return: protein stop gain classification for the list of matched tokens
        """
        gene_token, protein_stop_gain_token = tokens

        return ProteinStopGainClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos=protein_stop_gain_token.pos,
            ref=protein_stop_gain_token.ref,
            alt=protein_stop_gain_token.alt,
        )
