"""The module for Genomic Insertion Validation."""

from variation.schemas.classification_response_schema import (
    ClassificationType,
    GenomicInsertionClassification,
    Nomenclature,
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import GenomicValidator


class GenomicInsertion(GenomicValidator):
    """The Genomic Insertion Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicInsertionClassification, accessions: list[str]
    ) -> list[ValidationResult]:
        """Get list of validation results for a given classification and accessions

        :param classification: A classification for a list of tokens
        :param accessions: A list of accessions for a classification
        :return: List of validation results containing invalid and valid results
        """
        invalid_pos_msg = self.validate_5_prime_to_3_prime(
            classification.pos0, pos1=classification.pos1
        )
        if invalid_pos_msg:
            return [
                ValidationResult(
                    accession=None,
                    classification=classification,
                    is_valid=False,
                    errors=[invalid_pos_msg],
                )
            ]

        validation_results = []

        if classification.nomenclature == Nomenclature.GNOMAD_VCF:
            ref = classification.matching_tokens[0].ref
        else:
            ref = None

        for alt_ac in accessions:
            errors = []

            if ref:
                # gnomAD VCF provides reference, so we should validate this
                invalid_ref_msg = self.validate_reference_sequence(
                    alt_ac,
                    classification.pos0,
                    end_pos=classification.pos1,
                    expected_ref=ref,
                )
                if invalid_ref_msg:
                    errors.append(invalid_ref_msg)
            else:
                # Validate ac and pos
                invalid_ac_pos_msg = self.validate_ac_and_pos(
                    alt_ac, classification.pos0, end_pos=classification.pos1
                )
                if invalid_ac_pos_msg:
                    errors.append(invalid_ac_pos_msg)

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
        """Return whether or not the classification type is genomic insertion"""
        return classification_type == ClassificationType.GENOMIC_INSERTION
