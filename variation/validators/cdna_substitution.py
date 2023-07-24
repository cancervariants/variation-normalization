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
        self, classification: CdnaSubstitutionClassification, accessions: List[str]
    ) -> List[ValidationResult]:
        validation_results = []

        for c_ac in accessions:
            errors = []
            cds_start, cds_start_err_msg = await self.get_cds_start(c_ac)

            if cds_start_err_msg:
                errors.append(cds_start_err_msg)
            else:
                valid_ref_seq_msg = self.validate_reference_sequence(
                    c_ac, classification.pos + cds_start, classification.pos + cds_start,
                    classification.ref
                )
                if valid_ref_seq_msg:
                    errors.append(valid_ref_seq_msg)

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
        return "cdna substitution"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is cdna substitution."""
        return classification_type == ClassificationType.CDNA_SUBSTITUTION

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
