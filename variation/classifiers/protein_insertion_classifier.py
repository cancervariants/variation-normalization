"""A module for the Protein Insertion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, ProteinInsertionClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class ProteinInsertionClassifier(Classifier):
    """The Protein Insertion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Insertion classification type."""
        return ClassificationType.PROTEIN_INSERTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.PROTEIN_INSERTION]
        ]

    def match(self, tokens: List[Token]):
        gene_token, protein_ins_token = tokens

        return ProteinInsertionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            aa0=protein_ins_token.aa0,
            pos0=protein_ins_token.pos0,
            aa1=protein_ins_token.aa1,
            pos1=protein_ins_token.pos1,
            inserted_sequence=protein_ins_token.inserted_sequence
        )