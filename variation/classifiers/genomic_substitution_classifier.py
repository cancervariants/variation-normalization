"""A module for the Genomic Substitution Classifier."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, GenomicSubstitutionClassification, Nomenclature,
    SequenceOntology
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class GenomicSubstitutionClassifier(Classifier):
    """The Genomic Substitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Substitution classification type."""
        return ClassificationType.GENOMIC_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.GENOMIC_SUBSTITUTION]
        ]

    def match(self, tokens: List[Token]):
        gene_token, genomic_sub_token = tokens
        len_ref = len(genomic_sub_token.ref)
        len_alt = len(genomic_sub_token.alt)

        so_id = None
        if len_ref == 1 and len_alt == 1:
            so_id = SequenceOntology.SNV
        elif len_ref == len_alt:
            so_id = SequenceOntology.MNV

        if so_id:
            return GenomicSubstitutionClassification(
                matching_tokens=tokens,
                nomenclature=Nomenclature.FREE_TEXT,
                gene_token=gene_token,
                pos=genomic_sub_token.pos,
                ref=genomic_sub_token.ref,
                alt=genomic_sub_token.alt,
                so_id=so_id
            )
