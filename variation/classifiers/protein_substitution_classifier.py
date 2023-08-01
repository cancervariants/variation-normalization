"""A module for the Protein Substitution Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    Nomenclature,
    ProteinSubstitutionClassification,
)
from variation.schemas.token_response_schema import Token, TokenType


class ProteinSubstitutionClassifier(Classifier):
    """The ProteinSubstitution Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the protein substitution
        classification.

        :return: List of list of tokens, where order matters, that represent a protein
        substitution classification.
        """
        return [[TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION]]

    def match(self, tokens: List[Token]) -> ProteinSubstitutionClassification:
        """Return the protein substitution classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            protein substitution classification
        :return: protein substitution classification for the list of matched tokens
        """
        gene_token, protein_sub_token = tokens

        return ProteinSubstitutionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos=protein_sub_token.pos,
            ref=protein_sub_token.ref,
            alt=protein_sub_token.alt,
        )
