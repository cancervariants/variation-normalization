"""A module for the Cdna DelIns Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    CdnaDelInsClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType


class CdnaDelInsClassifier(Classifier):
    """The Cdna DelIns Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the cdna delins classification.

        :return: List of list of tokens, where order matters, that represent a cdna
        delins classification.
        """
        return [[TokenType.GENE, TokenType.CDNA_DELINS]]

    def match(self, tokens: List[Token]) -> CdnaDelInsClassification:
        """Return the cdna delins classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            cdna delins classification
        :return: cdna delins classification for the list of matched tokens
        """
        gene_token, cdna_delins_token = tokens

        return CdnaDelInsClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=cdna_delins_token.pos0,
            pos1=cdna_delins_token.pos1,
            inserted_sequence=cdna_delins_token.inserted_sequence,
        )
