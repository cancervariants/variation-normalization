"""A module for the DNA Coding Silent Mutation Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class CodingDNASilentMutationClassifier(SetBasedClassifier):
    """The Coding DNA Silent Mutation Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA Silent Mutation classification type."""
        return ClassificationType.CODING_DNA_SILENT_MUTATION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.CODING_DNA_SILENT_MUTATION],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.CODING_DNA_SILENT_MUTATION],  # noqa: E501
            [TokenType.CODING_DNA_SILENT_MUTATION, TokenType.GENE],
            [TokenType.GENE, TokenType.CODING_DNA_SILENT_MUTATION],
            [TokenType.HGVS, TokenType.CODING_DNA_SILENT_MUTATION],
            [TokenType.REFERENCE_SEQUENCE, TokenType.CODING_DNA_SILENT_MUTATION],
            [TokenType.LOCUS_REFERENCE_GENOMIC, TokenType.CODING_DNA_SILENT_MUTATION]
        ]
