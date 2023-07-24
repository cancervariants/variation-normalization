"""The module for Genomic DelIns Validation."""
from typing import List

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, GenomicDelInsClassification, Nomenclature
)
from variation.schemas.validation_response_schema import ValidationResult
from .validator import Validator


class GenomicDelIns(Validator):
    """The Genomic DelIns Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicDelInsClassification, accessions: List[str]
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
        # TODO: Validate pos0 and pos1 exist on given accession

        for alt_ac in accessions:
            validation_results.append(
                ValidationResult(
                    accession=alt_ac,
                    classification=classification,
                    is_valid=True,
                    errors=[]
                )
            )

        return validation_results

    def variation_name(self) -> str:
        """Return the variation name."""
        return "genomic delins"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is genomic delins"""
        return classification_type == ClassificationType.GENOMIC_DELINS

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
