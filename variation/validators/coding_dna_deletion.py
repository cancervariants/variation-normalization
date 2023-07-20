"""The module for cDNA Deletion Validation."""
from typing import Optional, List

from variation.schemas.classification_response_schema import (
    ClassificationType, Classification, CdnaDeletionClassification, Nomenclature
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class CdnaDeletion(Validator):
    """The cDNA Deletion Validator class."""

    async def get_valid_invalid_results(
        self, classification: CdnaDeletionClassification, transcripts: List[str]
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
            else:
                # validate deleted sequence
                if classification.nomenclature in {Nomenclature.FREE_TEXT,
                                                   Nomenclature.HGVS}:
                    # HGVS deleted sequence includes start and end
                    if classification.deleted_sequence:
                        invalid_del_seq_msg = self.validate_reference_sequence(
                            t, cds_start + classification.pos0,
                            cds_start + classification.pos1 + 1,
                            classification.deleted_sequence
                        )

                        if invalid_del_seq_msg:
                            errors.append(invalid_del_seq_msg)

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
        return "cdna deletion"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is cdna deletion."""
        return classification_type == ClassificationType.CODING_DNA_DELETION

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
