"""A module for the Genomic Substitution Classifier."""
from typing import List, Optional

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    GenomicSubstitutionClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType


class GenomicSubstitutionClassifier(Classifier):
    """The Genomic Substitution Classifier class."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the genomic substitution
        classification.

        :return: List of list of tokens, where order matters, that represent a genomic
        substitution classification.
        """
        return [
            [TokenType.GENE, TokenType.GENOMIC_SUBSTITUTION],
            [
                TokenType.GENE,
                TokenType.PROTEIN_SUBSTITUTION,
                TokenType.GENOMIC_SUBSTITUTION,
            ],
        ]

    def match(self, tokens: List[Token]) -> Optional[GenomicSubstitutionClassification]:
        """Return the genomic substitution classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            genomic substitution classification
        :return: genomic substitution classification for the list of matched tokens
        """
        if len(tokens) == 2:
            gene_token, genomic_sub_token = tokens
        else:
            gene_token, _, genomic_sub_token = tokens

        if len(genomic_sub_token.ref) == len(genomic_sub_token.alt):
            return GenomicSubstitutionClassification(
                matching_tokens=tokens,
                nomenclature=Nomenclature.FREE_TEXT,
                gene_token=gene_token,
                pos=genomic_sub_token.pos,
                ref=genomic_sub_token.ref,
                alt=genomic_sub_token.alt,
            )
