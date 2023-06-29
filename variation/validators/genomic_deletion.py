"""The module for Genomic Deletion Validation."""
from typing import List, Optional

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, GenomicDeletionClassification, Nomenclature
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from .validator import Validator


class GenomicDeletion(Validator):
    """The Genomic Deletion Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicDeletionClassification,
        transcripts: List[str], gene_tokens: List[GeneToken]
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
        # _validate_gene_pos?

        for ac in transcripts:
            errors = []

            invalid_ac_pos = self.validate_ac_and_pos(
                ac, classification.pos0, end_pos=classification.pos1
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
                            ac, classification.pos0, classification.pos1 + 1,
                            classification.deleted_sequence
                        )

                        if invalid_del_seq_message:
                            errors.append(invalid_del_seq_message)

            if not errors:
                if classification.nomenclature == Nomenclature.GNOMAD_VCF:
                    validate_ref_msg = self.validate_reference_sequence(
                        ac, classification.pos0, classification.pos1 + 1,
                        classification.matching_tokens[0].ref
                    )

                    if validate_ref_msg:
                        errors.append(validate_ref_msg)

            validation_results.append(
                ValidationResult(
                    accession=ac,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors,
                    gene_tokens=gene_tokens
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

    async def get_transcripts(
        self, gene_tokens: List, classification: Classification, errors: List
    ) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        if classification.nomenclature == Nomenclature.HGVS:
            transcripts = [classification.ac]
        else:
            transcripts = await self.get_genomic_transcripts(
                classification, gene_tokens, errors
            )
        return transcripts

    def get_gene_tokens(self, classification: Classification) -> List:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)
