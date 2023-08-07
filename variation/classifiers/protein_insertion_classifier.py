"""A module for the Protein Insertion Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    Nomenclature,
    ProteinInsertionClassification,
)
from variation.schemas.token_response_schema import Token, TokenType


class ProteinInsertionClassifier(Classifier):
    """The Protein Insertion Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the protein insertion classification.

        :return: List of list of tokens, where order matters, that represent a protein
        insertion classification.
        """
        return [[TokenType.GENE, TokenType.PROTEIN_INSERTION]]

    def match(self, tokens: List[Token]) -> ProteinInsertionClassification:
        """Return the protein insertion classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            protein insertion classification
        :return: protein insertion classification for the list of matched tokens
        """
        gene_token, protein_ins_token = tokens

        return ProteinInsertionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            aa0=protein_ins_token.aa0,
            pos0=protein_ins_token.pos0,
            aa1=protein_ins_token.aa1,
            pos1=protein_ins_token.pos1,
            inserted_sequence=protein_ins_token.inserted_sequence,
        )
