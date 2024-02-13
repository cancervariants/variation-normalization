"""A module for the Protein DelIns Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    Nomenclature,
    ProteinDelInsClassification,
)
from variation.schemas.token_response_schema import Token, TokenType


class ProteinDelInsClassifier(Classifier):
    """The Protein DelIns Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the protein delins classification.

        :return: List of list of tokens, where order matters, that represent a protein
        delins classification.
        """
        return [[TokenType.GENE, TokenType.PROTEIN_DELINS]]

    def match(self, tokens: List[Token]) -> ProteinDelInsClassification:
        """Return the protein delins classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            protein delins classification
        :return: protein delins classification for the list of matched tokens
        """
        gene_token, protein_delins_token = tokens

        return ProteinDelInsClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            aa0=protein_delins_token.aa0,
            pos0=protein_delins_token.pos0,
            aa1=protein_delins_token.aa1,
            pos1=protein_delins_token.pos1,
            inserted_sequence=protein_delins_token.inserted_sequence,
        )
