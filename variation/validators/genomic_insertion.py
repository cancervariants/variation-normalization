"""The module for Genomic Insertion Validation."""
from typing import List, Optional

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, GenomicInsertionClassification
)
from variation.schemas.token_response_schema import Token, TokenType, GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from .validator import Validator


class GenomicInsertion(Validator):
    """The Genomic Insertion Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicInsertionClassification,
        transcripts: List[str], gene_tokens: List[GeneToken]
    ) -> List[ValidationResult]:
        validation_results = []

        for t in transcripts:
            errors = []

            if (classification.pos0 >= classification.pos1):
                errors.append(
                    "Start position cannot be greater or equal than end position"
                )

            validation_results.append(
                ValidationResult(
                    accession=t,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors,
                    gene_tokens=gene_tokens
                )
            )

        return validation_results

    def variation_name(self) -> str:
        """Return the variation name."""
        return "genomic insertion"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is genomic insertion"""
        return classification_type == ClassificationType.GENOMIC_INSERTION

    async def get_transcripts(
        self, gene_tokens: List, classification: Classification, errors: List
    ) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        transcripts = await self.get_genomic_transcripts(
            classification, gene_tokens, errors
        )
        return transcripts

    def get_gene_tokens(self, classification: Classification) -> List:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)
