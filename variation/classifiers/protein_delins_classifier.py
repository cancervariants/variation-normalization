"""A module for the Protein DelIns Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, ProteinDelInsClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class ProteinDelInsClassifier(Classifier):
    """The Protein DelIns Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein DelIns classification type."""
        return ClassificationType.PROTEIN_DELINS

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.PROTEIN_DELINS]
        ]

    def match(self, tokens: List[Token]):
        gene_token, protein_delins_token = tokens

        return ProteinDelInsClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene=gene_token,
            aa0=protein_delins_token.aa0,
            pos0=protein_delins_token.pos0,
            aa1=protein_delins_token.aa1,
            pos1=protein_delins_token.pos1,
            inserted_sequence=protein_delins_token.inserted_sequence
        )
