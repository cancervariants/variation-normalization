"""The module for Genomic Reference Agree Validation."""
from typing import Optional, List

from variation.schemas.classification_response_schema import (
    ClassificationType, Classification, GenomicReferenceAgreeClassification,
    Nomenclature
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class GenomicReferenceAgree(Validator):
    """The Genomic Reference Agree Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicReferenceAgreeClassification,
        transcripts: List[str]
    ) -> List[ValidationResult]:
        validation_results = []

        for t in transcripts:
            errors = []

            # TODO: Validate pos exists on given accession
            if classification.nomenclature == Nomenclature.GNOMAD_VCF:
                token = classification.matching_tokens[0]
                ref = token.ref
                pos = token.pos

                invalid_ref_msg = self.validate_reference_sequence(t, pos, pos, ref)
                if invalid_ref_msg:
                    errors.append(invalid_ref_msg)

            validation_results.append(
                ValidationResult(
                    accession=t,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors
                )
            )

        return validation_results

    def variation_name(self) -> str:
        """Return the variation name."""
        return "genomic reference agree"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is genomic reference agree"""
        return classification_type == ClassificationType.GENOMIC_REFERENCE_AGREE

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
