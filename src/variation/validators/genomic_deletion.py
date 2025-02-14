"""The module for Genomic Deletion Validation."""

from variation.schemas.classification_response_schema import (
    Classification,
    ClassificationType,
    GenomicDeletionClassification,
    Nomenclature,
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class GenomicDeletion(Validator):
    """The Genomic Deletion Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicDeletionClassification, accessions: list[str]
    ) -> list[ValidationResult]:
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

        for alt_ac in accessions:
            errors = []

            invalid_ac_pos = self.validate_ac_and_pos(
                alt_ac, classification.pos0, end_pos=classification.pos1
            )
            if invalid_ac_pos:
                errors.append(invalid_ac_pos)
            else:
                if (
                    classification.nomenclature
                    in {
                        Nomenclature.FREE_TEXT,
                        Nomenclature.HGVS,
                    }
                    and classification.deleted_sequence
                ):
                    # Validate deleted sequence
                    # HGVS deleted sequence includes start and end
                    invalid_del_seq_message = self.validate_reference_sequence(
                        alt_ac,
                        classification.pos0,
                        classification.pos1
                        if classification.pos1
                        else classification.pos0,
                        classification.deleted_sequence,
                    )

                    if invalid_del_seq_message:
                        errors.append(invalid_del_seq_message)

            if not errors and classification.nomenclature == Nomenclature.GNOMAD_VCF:
                # Validate reference sequence
                ref = classification.matching_tokens[0].ref
                validate_ref_msg = self.validate_reference_sequence(
                    alt_ac,
                    classification.pos0 - 1,
                    end_pos=classification.pos0 + (len(ref) - 1),
                    expected_ref=ref,
                )

                if validate_ref_msg:
                    errors.append(validate_ref_msg)

            if not errors and classification.gene_token:
                # Validate positions exist within gene range
                invalid_gene_pos_msg = await self._validate_gene_pos(
                    classification.gene_token.matched_value,
                    alt_ac,
                    classification.pos0,
                    classification.pos1,
                )
                if invalid_gene_pos_msg:
                    errors.append(invalid_gene_pos_msg)

            validation_results.append(
                ValidationResult(
                    accession=alt_ac,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors,
                )
            )

        return validation_results

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is genomic deletion"""
        return classification_type == ClassificationType.GENOMIC_DELETION

    async def get_accessions(
        self, classification: Classification, errors: list
    ) -> list[str]:
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
            accessions = await self.get_genomic_accessions(classification, errors)
        return accessions
