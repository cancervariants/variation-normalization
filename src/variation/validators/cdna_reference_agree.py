"""The module for Cdna Substitution Validation."""

from variation.schemas.classification_response_schema import (
    CdnaReferenceAgreeClassification,
    Classification,
    ClassificationType,
    Nomenclature,
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class CdnaReferenceAgree(Validator):
    """The Cdna Reference Agree Validator class."""

    async def get_valid_invalid_results(
        self, classification: CdnaReferenceAgreeClassification, accessions: list[str]
    ) -> list[ValidationResult]:
        """Get list of validation results for a given classification and accessions

        :param classification: A classification for a list of tokens
        :param accessions: A list of accessions for a classification
        :return: List of validation results containing invalid and valid results
        """
        validation_results = []

        for c_ac in accessions:
            errors = []
            cds_start, cds_start_err_msg = await self.get_cds_start(c_ac)

            if cds_start_err_msg:
                errors.append(cds_start_err_msg)
            else:
                invalid_ac_pos_msg = self.validate_ac_and_pos(
                    c_ac, cds_start + classification.pos
                )
                if invalid_ac_pos_msg:
                    errors.append(invalid_ac_pos_msg)

            validation_results.append(
                ValidationResult(
                    accession=c_ac,
                    classification=classification,
                    cds_start=cds_start,
                    is_valid=not errors,
                    errors=errors,
                )
            )

        return validation_results

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is cdna reference agree."""
        return classification_type == ClassificationType.CDNA_REFERENCE_AGREE

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
            accessions = self.get_cdna_accessions(classification.gene_token, errors)
        return accessions
