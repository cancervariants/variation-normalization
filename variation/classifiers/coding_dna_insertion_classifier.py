"""A module for the Coding DNA insertion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, CdnaInsertionClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class CodingDNAInsertionClassifier(Classifier):
    """The Coding DNA insertion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA insertion classification type."""
        return ClassificationType.CODING_DNA_INSERTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.CODING_DNA_INSERTION]
        ]

    def match(self, tokens: List[Token]) -> CdnaInsertionClassification:
        gene_token, cdna_ins_token = tokens

        return CdnaInsertionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene=gene_token,
            pos0=cdna_ins_token.pos0,
            pos1=cdna_ins_token.pos1,
            inserted_sequence=cdna_ins_token.inserted_sequence
        )
