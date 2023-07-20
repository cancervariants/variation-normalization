"""The module for Coding DNA Insertion Validation."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, Classification, CdnaInsertionClassification, Nomenclature
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class CdnaInsertion(Validator):
    """The Coding DNA Insertion Validator class."""

    async def get_valid_invalid_results(
        self, classification: CdnaInsertionClassification, transcripts: List[str]
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

        for t in transcripts:
            errors = []
            cds_start, cds_start_err_msg = await self.get_cds_start(t)

            if cds_start_err_msg:
                errors.append(cds_start_err_msg)

            # TODO: Validate pos0 and pos1 exist on given accession

            validation_results.append(
                ValidationResult(
                    accession=t,
                    classification=classification,
                    cds_start=cds_start,
                    is_valid=not errors,
                    errors=errors
                )
            )

        return validation_results

    def variation_name(self) -> str:
        """Return the variation name."""
        return "cdna insertion"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is cdna insertion."""
        return classification_type == ClassificationType.CODING_DNA_INSERTION

    async def get_transcripts(
        self, classification: Classification, errors: List
    ) -> List[str]:
        """Get transcript accessions for a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        if classification.nomenclature == Nomenclature.HGVS:
            transcripts = [classification.ac]
        else:
            transcripts = self.get_coding_dna_transcripts(
                classification.gene_token, errors
            )
        return transcripts
