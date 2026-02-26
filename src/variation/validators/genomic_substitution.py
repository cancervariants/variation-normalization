"""The module for Genomic Substitution Validation."""

from variation.schemas.classification_response_schema import (
    ClassificationType,
    GenomicSubstitutionClassification,
    Nomenclature,
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import GenomicValidator


class GenomicSubstitution(GenomicValidator):
    """The Genomic Substitution Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicSubstitutionClassification, accessions: list[str]
    ) -> list[ValidationResult]:
        """Get list of validation results for a given classification and accessions

        :param classification: A classification for a list of tokens
        :param accessions: A list of accessions for a classification
        :return: List of validation results containing invalid and valid results
        """
        validation_results = []

        if classification.nomenclature == Nomenclature.GNOMAD_VCF:
            end_pos = classification.pos + (len(classification.alt) - 1)
        else:
            # HGVS is only 1 nuc
            end_pos = classification.pos

        for alt_ac in accessions:
            errors = []

            valid_ref_seq_msg = self.validate_reference_sequence(
                alt_ac, classification.pos, end_pos, classification.ref
            )
            if valid_ref_seq_msg:
                errors.append(valid_ref_seq_msg)

            validation_results.append(
                ValidationResult(
                    accession=alt_ac,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors,
                )
            )

        return validation_results

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is genomic
        substitution.
        """
        return classification_type == ClassificationType.GENOMIC_SUBSTITUTION
