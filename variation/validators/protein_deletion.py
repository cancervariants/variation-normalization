"""The module for Protein Deletion Validation."""
from typing import List, Optional

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, Nomenclature, ProteinDeletionClassification
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator
from variation.utils import get_aa1_codes


class ProteinDeletion(Validator):
    """The Protein Deletion Validator class."""

    async def get_valid_invalid_results(
        self, classification: ProteinDeletionClassification, transcripts: List[str]
    ) -> List[ValidationResult]:
        errors = []
        if classification.pos1:
            if classification.pos0 >= classification.pos1:
                errors.append(
                    "Positions deleted should contain two different positions and "
                    "should be listed from 5' to 3'"
                )

        # Only HGVS Expressions are validated
        # Free text is validated during tokenization
        if classification.nomenclature == Nomenclature.HGVS:
            # 1 letter AA codes for aa0
            aa0_codes = get_aa1_codes(classification.aa0)
            if aa0_codes:
                classification.aa0 = aa0_codes
            else:
                errors.append(f"`aa0` not valid amino acid(s): {classification.aa0}")

            if classification.aa1:
                # 1 letter AA codes for aa1
                aa1_codes = get_aa1_codes(classification.aa1)
                if aa1_codes:
                    classification.aa1 = aa1_codes
                else:
                    errors.append(
                        f"`aa1` not valid amino acid(s): {classification.aa1}"
                    )

        if errors:
            return [
                ValidationResult(
                    accession=None,
                    classification=classification,
                    is_valid=False,
                    errors=errors
                )
            ]

        validation_results = []

        for p_ac in transcripts:
            errors = []

            # Validate aa0 exists at pos0 on given protein accession
            invalid_aa0_seq_msg = self.validate_reference_sequence(
                p_ac, classification.pos0, classification.pos0, classification.aa0
            )
            if invalid_aa0_seq_msg:
                errors.append(invalid_aa0_seq_msg)

            # Validate aa1 exists at pos1
            if classification.aa1 and classification.pos1:
                invalid_aa1_seq_msg = self.validate_reference_sequence(
                    p_ac, classification.pos1, classification.pos1, classification.aa1
                )

                if invalid_aa1_seq_msg:
                    errors.append(invalid_aa1_seq_msg)

            # Validate that deleted sequence matches expected
            if classification.nomenclature in {Nomenclature.FREE_TEXT,
                                               Nomenclature.HGVS}:
                # HGVS deleted sequence includes start and end
                if classification.deleted_sequence:
                    if classification.pos1 is not None:
                        invalid_del_seq_msg = self.validate_reference_sequence(
                            p_ac, classification.pos0, classification.pos1 + 1,
                            classification.deleted_sequence
                        )

                        if invalid_del_seq_msg:
                            errors.append(invalid_del_seq_msg)

            validation_results.append(
                ValidationResult(
                    accession=p_ac,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors
                )
            )

        return validation_results

    def variation_name(self) -> str:
        """Return the variation name."""
        return "protein deletion"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is protein deletion."""
        return classification_type == ClassificationType.PROTEIN_DELETION

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
            transcripts = self.get_protein_transcripts(
                classification.gene_token, errors
            )
        return transcripts
