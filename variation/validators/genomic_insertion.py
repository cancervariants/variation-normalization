"""The module for Genomic Insertion Validation."""
from typing import List

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, GenomicInsertionClassification, Nomenclature
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class GenomicInsertion(Validator):
    """The Genomic Insertion Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicInsertionClassification, transcripts: List[str]
    ) -> List[ValidationResult]:
        if classification.pos1 and classification.pos0 >= classification.pos1:
            return [ValidationResult(
                accession=None,
                classification=classification,
                is_valid=False,
                errors=[(
                    "Positions deleted should contain two different positions and "
                    "should be listed from 5' to 3'")]
            )]

        validation_results = []

        for ac in transcripts:
            errors = []

            # gnomAD VCF provides reference, so we should validate this
            if classification.nomenclature == Nomenclature.GNOMAD_VCF:
                invalid_ref_msg = self.validate_reference_sequence(
                    ac, classification.pos0, classification.pos1,
                    classification.matching_tokens[0].ref
                )
                if invalid_ref_msg:
                    errors.append(invalid_ref_msg)

            validation_results.append(
                ValidationResult(
                    accession=ac,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors
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
        self, classification: Classification, errors: List
    ) -> List[str]:
        """Get transcript accessions for a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        if classification.nomenclature == Nomenclature.HGVS:
            transcripts = [classification.ac]
        else:
            transcripts = await self.get_genomic_transcripts(
                classification, errors
            )
        return transcripts
