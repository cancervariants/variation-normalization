"""The module for Genomic Duplication Validation."""
from typing import List, Optional

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, GenomicDuplicationClassification, Nomenclature
)
from variation.schemas.validation_response_schema import ValidationResult
from .validator import Validator


class GenomicDuplication(Validator):
    """The Genomic Duplication Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicDuplicationClassification, accessions: List[str]
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
            errors = []

            if classification.gene_token:
                invalid_gene_pos_msg = await self._validate_gene_pos(
                    classification.gene_token.matched_value, alt_ac,
                    classification.pos0, classification.pos1
                )
                if invalid_gene_pos_msg:
                    errors.append(invalid_gene_pos_msg)

            validation_results.append(
                ValidationResult(
                    accession=alt_ac,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors
                )
            )

        return validation_results

    def variation_name(self) -> str:
        """Return the variation name."""
        return "genomic duplication"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is genomic duplication"""
        return classification_type == ClassificationType.GENOMIC_DUPLICATION

    async def get_accessions(
        self, classification: Classification, errors: List
    ) -> Optional[List[str]]:
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
