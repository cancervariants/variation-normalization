"""The module for Genomic Duplication Validation."""
from typing import List, Optional

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, GenomicDuplicationClassification, Nomenclature
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from .validator import Validator


class GenomicDuplication(Validator):
    """The Genomic Duplication Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicDuplicationClassification, transcripts: List[str]
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
        # _validate_gene_pos?

        for ac in transcripts:
            validation_results.append(
                ValidationResult(
                    accession=ac,
                    classification=classification,
                    is_valid=True,
                    errors=[]
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

    async def get_transcripts(
        self, classification: Classification, errors: List
    ) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        if classification.nomenclature == Nomenclature.HGVS:
            transcripts = [classification.ac]
        else:
            transcripts = await self.get_genomic_transcripts(
                classification, errors
            )
        return transcripts
