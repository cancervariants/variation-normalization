"""The module for Genomic Deletion Validation."""
from typing import List

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, GenomicDeletionClassification, Nomenclature
)
from variation.schemas.validation_response_schema import ValidationResult
from .validator import Validator


class GenomicDeletion(Validator):
    """The Genomic Deletion Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicDeletionClassification, accessions: List[str]
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

            invalid_ac_pos = self.validate_ac_and_pos(
                alt_ac, classification.pos0, end_pos=classification.pos1
            )
            if invalid_ac_pos:
                errors.append(invalid_ac_pos)
            else:
                # validate deleted sequence
                if classification.nomenclature in {Nomenclature.FREE_TEXT,
                                                   Nomenclature.HGVS}:
                    # HGVS deleted sequence includes start and end
                    if classification.deleted_sequence:
                        invalid_del_seq_message = self.validate_reference_sequence(
                            alt_ac, classification.pos0, classification.pos1 + 1,
                            classification.deleted_sequence
                        )

                        if invalid_del_seq_message:
                            errors.append(invalid_del_seq_message)

            if not errors:
                if classification.nomenclature == Nomenclature.GNOMAD_VCF:
                    validate_ref_msg = self.validate_reference_sequence(
                        alt_ac, classification.pos0 - 1, classification.pos1,
                        classification.matching_tokens[0].ref
                    )

                    if validate_ref_msg:
                        errors.append(validate_ref_msg)

            if not errors and classification.gene_token:
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
        return "genomic deletion"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is genomic deletion"""
        return classification_type == ClassificationType.GENOMIC_DELETION

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
