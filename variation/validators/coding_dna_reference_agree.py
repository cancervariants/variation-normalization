"""The module for Coding DNA Substitution Validation."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, Classification, Nomenclature, CdnaReferenceAgreeClassification
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class CdnaReferenceAgree(Validator):
    """The Coding DNA Reference Agree Validator class."""

    async def get_valid_invalid_results(
        self, classification: CdnaReferenceAgreeClassification, accessions: List[str]
    ) -> List[ValidationResult]:
        validation_results = []

        for c_ac in accessions:
            errors = []
            cds_start, cds_start_err_msg = await self.get_cds_start(c_ac)

            if cds_start_err_msg:
                errors.append(cds_start_err_msg)

            # TODO: Validate pos exists on given accession

            validation_results.append(
                ValidationResult(
                    accession=c_ac,
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

    async def get_accessions(
        self, classification: Classification, errors: List
    ) -> List[str]:
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
            accessions = self.get_cdna_accessions(
                classification.gene_token, errors
            )
        return accessions
