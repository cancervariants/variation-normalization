"""A module for the Cdna Substitution Classifier."""
from typing import List, Optional

from variation.schemas.classification_response_schema import (
    ClassificationType, CdnaSubstitutionClassification, Nomenclature, SequenceOntology
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class CdnaSubstitutionClassifier(Classifier):
    """The Cdna Substitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Cdna Substitution classification type."""
        return ClassificationType.CDNA_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.CDNA_SUBSTITUTION],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.CDNA_SUBSTITUTION]  # noqa: E501
        ]

    def match(self, tokens: List[Token]) -> Optional[CdnaSubstitutionClassification]:
        if len(tokens) == 2:
            gene_token, cdna_sub_token = tokens
        else:
            gene_token, _, cdna_sub_token = tokens

        len_ref = len(cdna_sub_token.ref)
        len_alt = len(cdna_sub_token.alt)

        so_id = None
        if len_ref == 1 and len_alt == 1:
            so_id = SequenceOntology.SNV
        elif len_ref == len_alt:
            so_id = SequenceOntology.MNV

        if so_id:
            return CdnaSubstitutionClassification(
                matching_tokens=tokens,
                nomenclature=Nomenclature.FREE_TEXT,
                gene_token=gene_token,
                pos=cdna_sub_token.pos,
                ref=cdna_sub_token.ref,
                alt=cdna_sub_token.alt,
                so_id=so_id
            )
