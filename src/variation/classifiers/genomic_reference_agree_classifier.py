"""A module for the Genomic Reference Agree Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    GenomicReferenceAgreeClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType


class GenomicReferenceAgreeClassifier(Classifier):
    """The Genomic Reference Agree Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the genomic reference agree
        classification.

        :return: List of list of tokens, where order matters, that represent a genomic
        reference agree classification.
        """
        return [[TokenType.GENE, TokenType.GENOMIC_REFERENCE_AGREE]]

    def match(self, tokens: List[Token]) -> GenomicReferenceAgreeClassification:
        """Return the genomic reference agree classification from a list of token
        matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            genomic reference agree classification
        :return: genomic reference agree classification for the list of matched tokens
        """
        gene_token, genomic_ref_agree_token = tokens

        return GenomicReferenceAgreeClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos=genomic_ref_agree_token.pos,
        )
