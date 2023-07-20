"""The module for Coding DNA Substitution Validation."""
from typing import List, Optional

from variation.schemas.classification_response_schema import (
    ClassificationType, Classification, Nomenclature, CdnaReferenceAgreeClassification
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class CdnaReferenceAgree(Validator):
    """The Coding DNA Reference Agree Validator class."""

    async def get_valid_invalid_results(
        self, classification: CdnaReferenceAgreeClassification, transcripts: List[str]
    ) -> List[ValidationResult]:
        validation_results = []

        for t in transcripts:
            errors = []
            cds_start, cds_start_err_msg = await self.get_cds_start(t)

            if cds_start_err_msg:
                errors.append(cds_start_err_msg)

            # TODO: Validate pos exists on given accession

            validation_results.append(
                ValidationResult(
                    accession=t,
                    classification=classification,
                    cds_start=cds_start,
                    is_valid=not errors,
                    errors=errors
                )
            )

        return validation_results

    def variation_name(self) -> str:
        """Return the variation name."""
        return "cdna reference agree"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is cdna reference agree."""
        return classification_type == ClassificationType.CODING_DNA_REFERENCE_AGREE

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
            transcripts = self.get_coding_dna_transcripts(
                classification.gene_token, errors
            )
        return transcripts
