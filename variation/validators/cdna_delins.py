"""The module for Cdna DelIns Validation."""
from typing import List, Optional, Union

from variation.schemas.classification_response_schema import (
    CdnaDelInsClassification,
    Classification,
    ClassificationType,
    Nomenclature,
)
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class CdnaDelIns(Validator):
    """The Cdna DelIns Validator class."""

    async def get_valid_invalid_results(
        self, classification: CdnaDelInsClassification, accessions: List[str]
    ) -> List[ValidationResult]:
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

        for c_ac in accessions:
            errors = []
            cds_start, cds_start_err_msg = await self.get_cds_start(c_ac)

            if cds_start_err_msg:
                errors.append(cds_start_err_msg)
            else:
                # Validate accession and positions
                invalid_ac_pos_msg = self.validate_ac_and_pos(
                    c_ac,
                    cds_start + classification.pos0,
                    end_pos=cds_start + classification.pos1
                    if classification.pos1
                    else None,
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
        """Return whether or not the classification type is cdna delins."""
        return classification_type == ClassificationType.CDNA_DELINS

    async def get_accessions(
        self,
        classification: Classification,
        errors: List,
        input_assembly: Optional[
            Union[ClinVarAssembly.GRCH37, ClinVarAssembly.GRCH38]
        ] = None,
    ) -> List[str]:
        """Get accessions for a given classification.
        If `classification.nomenclature == Nomenclature.HGVS`, will return the accession
        in the HGVS expression.
        Else, will get all accessions associated to the gene

        :param classification: The classification for list of tokens
        :param errors: List of errors
        :param input_assembly: Assembly used for initial input query. Only used when
            initial query is using genomic free text or gnomad vcf format
        :return: List of accessions
        """
        if classification.nomenclature == Nomenclature.HGVS:
            accessions = [classification.ac]
        else:
            accessions = self.get_cdna_accessions(classification.gene_token, errors)
        return accessions
