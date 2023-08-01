"""A module for the cDNA Reference Agree Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    CdnaReferenceAgreeClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType


class CdnaReferenceAgreeClassifier(Classifier):
    """The Cdna Reference Agree Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the cdna reference agree
        classification.

        :return: List of list of tokens, where order matters, that represent a cdna
        reference agree classification.
        """
        return [
            [TokenType.GENE, TokenType.CDNA_REFERENCE_AGREE],
            [
                TokenType.GENE,
                TokenType.PROTEIN_SUBSTITUTION,
                TokenType.CDNA_REFERENCE_AGREE,
            ],
        ]

    def match(self, tokens: List[Token]) -> CdnaReferenceAgreeClassification:
        """Return the cdna reference agree classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            cdna reference agree classification
        :return: cdna reference agree classification for the list of matched tokens
        """
        if len(tokens) == 2:
            gene_token, cdna_ref_agree_token = tokens
        else:
            gene_token, _, cdna_ref_agree_token = tokens

        return CdnaReferenceAgreeClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos=cdna_ref_agree_token.pos,
        )
