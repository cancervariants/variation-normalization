"""A module for the Genomic Silent Mutation Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class GenomicSilentMutationClassifier(SetBasedClassifier):
    """The Genomic Silent Mutation Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Silent Mutation classification type."""
        return ClassificationType.GENOMIC_SILENT_MUTATION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.CHROMOSOME, TokenType.GENOMIC_SILENT_MUTATION],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.GENOMIC_SILENT_MUTATION],  # noqa: E501
            [TokenType.GENOMIC_SILENT_MUTATION, TokenType.GENE],
            [TokenType.GENE, TokenType.GENOMIC_SILENT_MUTATION],
            [TokenType.HGVS, TokenType.GENOMIC_SILENT_MUTATION],
            [TokenType.REFERENCE_SEQUENCE, TokenType.GENOMIC_SILENT_MUTATION]
        ]
