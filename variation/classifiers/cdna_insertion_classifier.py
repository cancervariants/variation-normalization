"""A module for the Cdna insertion Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    CdnaInsertionClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType


class CdnaInsertionClassifier(Classifier):
    """The Cdna insertion Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the cdna insertion classification.

        :return: List of list of tokens, where order matters, that represent a cdna
        insertion classification.
        """
        return [[TokenType.GENE, TokenType.CDNA_INSERTION]]

    def match(self, tokens: List[Token]) -> CdnaInsertionClassification:
        """Return the cdna insertion classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            cdna insertion classification
        :return: cdna insertion classification for the list of matched tokens
        """
        gene_token, cdna_ins_token = tokens

        return CdnaInsertionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=cdna_ins_token.pos0,
            pos1=cdna_ins_token.pos1,
            inserted_sequence=cdna_ins_token.inserted_sequence,
        )
