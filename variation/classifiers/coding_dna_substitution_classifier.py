"""A module for the Coding DNA Substitution Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, CdnaSubstitutionClassification, Nomenclature, SequenceOntology
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class CodingDNASubstitutionClassifier(Classifier):
    """The Coding DNA Substitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA Substitution classification type."""
        return ClassificationType.CODING_DNA_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.CODING_DNA_SUBSTITUTION]
        ]

    def match(self, tokens: List[Token]):
        gene_token, cdna_sub_token = tokens
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
                gene=gene_token,
                pos=cdna_sub_token.pos,
                ref=cdna_sub_token.ref,
                alt=cdna_sub_token.alt,
                so_id=so_id
            )
