"""The module for cDNA Substitution Validation."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, Classification, CdnaSubstitutionClassification, Nomenclature
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class CdnaSubstitution(Validator):
    """The cDNA Substitution Validator class."""

    async def get_valid_invalid_results(
        self, classification: CdnaSubstitutionClassification, transcripts: List[str]
    ) -> List[ValidationResult]:
        validation_results = []

        for t in transcripts:
            errors = []
            cds_start, cds_start_err_msg = await self.get_cds_start(t)

            if cds_start_err_msg:
                errors.append(cds_start_err_msg)
            else:
                valid_ref_seq_msg = self.validate_reference_sequence(
                    t, classification.pos + cds_start, classification.pos + cds_start,
                    classification.ref
                )
                if valid_ref_seq_msg:
                    errors.append(valid_ref_seq_msg)

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
        return "cdna substitution"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is cdna substitution."""
        return classification_type == ClassificationType.CODING_DNA_SUBSTITUTION

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
