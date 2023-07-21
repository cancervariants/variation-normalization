"""A module for the Protein Deletion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, ProteinDeletionClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class ProteinDeletionClassifier(Classifier):
    """The Protein Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Deletion classification type."""
        return ClassificationType.PROTEIN_DELETION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.PROTEIN_DELETION]
        ]

    def match(self, tokens: List[Token]) -> ProteinDeletionClassification:
        gene_token, protein_del_token = tokens

        return ProteinDeletionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            aa0=protein_del_token.aa0,
            pos0=protein_del_token.pos0,
            aa1=protein_del_token.aa1,
            pos1=protein_del_token.pos1,
            deleted_sequence=protein_del_token.deleted_sequence
        )
