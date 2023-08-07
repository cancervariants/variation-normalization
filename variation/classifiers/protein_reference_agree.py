"""A module for the Reference Agree Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    Nomenclature,
    ProteinReferenceAgreeClassification,
)
from variation.schemas.token_response_schema import Token, TokenType


class ProteinReferenceAgreeClassifier(Classifier):
    """The Reference Agree Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the protein reference agree
        classification.

        :return: List of list of tokens, where order matters, that represent a protein
        reference agree classification.
        """
        return [[TokenType.GENE, TokenType.PROTEIN_REFERENCE_AGREE]]

    def match(self, tokens: List[Token]) -> ProteinReferenceAgreeClassification:
        """Return the protein reference agree classification from a list of token
        matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            protein reference agree classification
        :return: protein reference agree classification for the list of matched tokens
        """
        gene_token, protein_ref_agree_token = tokens

        return ProteinReferenceAgreeClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos=protein_ref_agree_token.pos,
            ref=protein_ref_agree_token.ref,
        )
