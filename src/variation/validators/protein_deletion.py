"""The module for Protein Deletion Validation."""

from variation.schemas.classification_response_schema import (
    Classification,
    ClassificationType,
    Nomenclature,
    ProteinDeletionClassification,
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class ProteinDeletion(Validator):
    """The Protein Deletion Validator class."""

    async def get_valid_invalid_results(
        self, classification: ProteinDeletionClassification, accessions: list[str]
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

        # Only HGVS Expressions are validated
        # Free text is validated during tokenization
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

            # Validate aa0 exists at pos0 on given protein accession
            invalid_aa0_seq_msg = self.validate_reference_sequence(
                p_ac, classification.pos0, classification.pos0, classification.aa0
            )
            if invalid_aa0_seq_msg:
                errors.append(invalid_aa0_seq_msg)

            # Validate aa1 exists at pos1
            if classification.aa1 and classification.pos1:
                invalid_aa1_seq_msg = self.validate_reference_sequence(
                    p_ac, classification.pos1, classification.pos1, classification.aa1
                )

                if invalid_aa1_seq_msg:
                    errors.append(invalid_aa1_seq_msg)

            # Validate that deleted sequence matches expected
            if (
                classification.nomenclature
                in {
                    Nomenclature.FREE_TEXT,
                    Nomenclature.HGVS,
                }
                and classification.deleted_sequence
                and classification.pos1 is not None
            ):
                # HGVS deleted sequence includes start and end
                invalid_del_seq_msg = self.validate_reference_sequence(
                    p_ac,
                    classification.pos0,
                    classification.pos1,
                    classification.deleted_sequence,
                )

                if invalid_del_seq_msg:
                    errors.append(invalid_del_seq_msg)

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
        """Return whether or not the classification type is protein deletion."""
        return classification_type == ClassificationType.PROTEIN_DELETION

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
