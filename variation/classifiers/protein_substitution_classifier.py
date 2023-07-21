"""A module for the Protein Substitution Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, ProteinSubstitutionClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class ProteinSubstitutionClassifier(Classifier):
    """The ProteinSubstitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Substitution classification type."""
        return ClassificationType.PROTEIN_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION]
        ]

    def match(self, tokens: List[Token]) -> ProteinSubstitutionClassification:
        gene_token, protein_sub_token = tokens

        return ProteinSubstitutionClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos=protein_sub_token.pos,
            ref=protein_sub_token.ref,
            alt=protein_sub_token.alt
        )
