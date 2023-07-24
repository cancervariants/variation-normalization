"""The module for Coding DNA DelIns Validation."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, Classification, Nomenclature, CdnaDelInsClassification
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class CodingDNADelIns(Validator):
    """The Coding DNA DelIns Validator class."""

    async def get_valid_invalid_results(
        self, classification: CdnaDelInsClassification, accessions: List[str]
    ) -> List[ValidationResult]:
        if classification.pos1 and classification.pos0 >= classification.pos1:
            return [ValidationResult(
                accession=None,
                classification=classification,
                is_valid=False,
                errors=[(
                    "Positions deleted should contain two different positions and "
                    "should be listed from 5' to 3'")]
            )]

        validation_results = []

        for c_ac in accessions:
            errors = []
            cds_start, cds_start_err_msg = await self.get_cds_start(c_ac)

            if cds_start_err_msg:
                errors.append(cds_start_err_msg)

            # TODO: Validate pos0 and pos1 exist on given accession

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
        return "cdna delins"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is cdna delins."""
        return classification_type == ClassificationType.CODING_DNA_DELINS

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
