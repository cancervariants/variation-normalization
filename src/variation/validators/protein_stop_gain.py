"""The module for Protein Stop Gain Validation."""

from variation.schemas.classification_response_schema import (
    Classification,
    ClassificationType,
    Nomenclature,
    ProteinStopGainClassification,
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class ProteinStopGain(Validator):
    """The Protein Stop Gain Validator class."""

    async def get_valid_invalid_results(
        self, classification: ProteinStopGainClassification, accessions: list[str]
    ) -> list[ValidationResult]:
        """Get list of validation results for a given classification and accessions

        :param classification: A classification for a list of tokens
        :param accessions: A list of accessions for a classification
        :return: List of validation results containing invalid and valid results
        """
        errors = []

        # Only HGVS Expressions are validated
        # Free text is validated during tokenization
        # Don't need to validate alt, since we know it's '*'
        if classification.nomenclature == Nomenclature.HGVS:
            invalid_classification_msgs = self.validate_protein_hgvs_classification(
                classification
            )
            if invalid_classification_msgs:
                return [
                    ValidationResult(
                        accession=None,
                        classification=classification,
                        is_valid=False,
                        errors=invalid_classification_msgs,
                    )
                ]

        validation_results = []

        for p_ac in accessions:
            errors = []

            valid_ref_seq_msg = self.validate_reference_sequence(
                p_ac, classification.pos, classification.pos, classification.ref
            )
            if valid_ref_seq_msg:
                errors.append(valid_ref_seq_msg)

            validation_results.append(
                ValidationResult(
                    accession=p_ac,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors,
                )
            )

        return validation_results

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is protein stop gain."""
        return classification_type == ClassificationType.PROTEIN_STOP_GAIN

    async def get_accessions(
        self, classification: Classification, errors: list
    ) -> list[str]:
        """Get accessions for a given classification.
        If `classification.nomenclature == Nomenclature.HGVS`, will return the accession
        in the HGVS expression.
        Else, will get all accessions associated to the gene

        :param classification: The classification for list of tokens
        :param errors: List of errors
        :return: List of accessions
        """
        if classification.nomenclature == Nomenclature.HGVS:
            accessions = [classification.ac]
        else:
            accessions = self.get_protein_accessions(classification.gene_token, errors)
        return accessions
