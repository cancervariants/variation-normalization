"""A module for the Protein Stop Gain Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    Nomenclature,
    ProteinStopGainClassification,
)
from variation.schemas.token_response_schema import Token, TokenType


class ProteinStopGainClassifier(Classifier):
    """The Protein Stop Gain Classifier class."""

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
