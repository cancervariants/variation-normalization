"""A module for the Genomic Duplication Ambiguous Classifier."""
from typing import List

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    GenomicDuplicationAmbiguousClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.utils import get_ambiguous_type


class GenomicDuplicationAmbiguousClassifier(Classifier):
    """The Genomic Duplication Ambiguous Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the genomic ambiguous duplication
        classification.

        :return: List of list of tokens, where order matters, that represent a genomic
        ambiguous duplication classification.
        """
        return [[TokenType.GENE, TokenType.GENOMIC_DUPLICATION_AMBIGUOUS]]

    def match(self, tokens: List[Token]) -> GenomicDuplicationAmbiguousClassification:
        """Return the genomic ambiguous duplication classification from a list of token
        matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            genomic ambiguous duplication classification
        :return: genomic ambiguous duplication classification for the list of matched
            tokens
        """
        gene_token, genomic_dup_token = tokens
        ambiguous_type = get_ambiguous_type(
            genomic_dup_token.pos0,
            genomic_dup_token.pos1,
            genomic_dup_token.pos2,
            genomic_dup_token.pos3,
            genomic_dup_token.ambiguous_regex_type,
        )

        return GenomicDuplicationAmbiguousClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=genomic_dup_token.pos0,
            pos1=genomic_dup_token.pos1,
            pos2=genomic_dup_token.pos2,
            pos3=genomic_dup_token.pos3,
            ambiguous_type=ambiguous_type,
        )
