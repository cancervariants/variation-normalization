"""The module for Genomic Insertion Validation."""
from typing import List

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, GenomicInsertionClassification, Nomenclature
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class GenomicInsertion(Validator):
    """The Genomic Insertion Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicInsertionClassification, accessions: List[str]
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

        for alt_ac in accessions:
            errors = []

            if classification.nomenclature == Nomenclature.GNOMAD_VCF:
                # gnomAD VCF provides reference, so we should validate this
                invalid_ref_msg = self.validate_reference_sequence(
                    alt_ac, classification.pos0, classification.pos1,
                    classification.matching_tokens[0].ref
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
                    errors=errors
                )
            )

        return validation_results

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is genomic insertion"""
        return classification_type == ClassificationType.GENOMIC_INSERTION

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
            accessions = await self.get_genomic_accessions(
                classification, errors
            )
        return accessions
