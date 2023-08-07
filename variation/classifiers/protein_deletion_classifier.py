"""A module for the Protein Deletion Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    Nomenclature,
    ProteinDeletionClassification,
)
from variation.schemas.token_response_schema import Token, TokenType


class ProteinDeletionClassifier(Classifier):
    """The Protein Deletion Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the protein deletion classification.

        :return: List of list of tokens, where order matters, that represent a protein
        deletion classification.
        """
        return [[TokenType.GENE, TokenType.PROTEIN_DELETION]]

    def match(self, tokens: List[Token]) -> ProteinDeletionClassification:
        """Return the protein deletion classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            protein deletion classification
        :return: protein deletion classification for the list of matched tokens
        """
        gene_token, protein_del_token = tokens

        return ProteinDeletionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            aa0=protein_del_token.aa0,
            pos0=protein_del_token.pos0,
            aa1=protein_del_token.aa1,
            pos1=protein_del_token.pos1,
            deleted_sequence=protein_del_token.deleted_sequence,
        )
