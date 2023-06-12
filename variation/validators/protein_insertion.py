"""The module for Protein Insertion Validation."""
from typing import List, Optional

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, Nomenclature, ProteinInsertionClassification
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator
from variation.utils import get_aa1_codes


class ProteinInsertion(Validator):
    """The Protein Insertion Validator class."""

    async def get_valid_invalid_results(
        self, classification: ProteinInsertionClassification,
        transcripts: List[str], gene_tokens: List[GeneToken]
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

            if classification.inserted_sequence:
                # 1 letter AA codes for inserted sequence
                ins_codes = get_aa1_codes(classification.inserted_sequence)
                if ins_codes:
                    classification.inserted_sequence = ins_codes
                else:
                    errors.append(
                        f"`inserted_sequence` not valid amino acid(s): {classification.inserted_sequence}"
                    )

        if errors:
            return [
                ValidationResult(
                    accession=None,
                    classification=classification,
                    is_valid=False,
                    errors=errors,
                    gene_tokens=gene_tokens
                )
            ]

        validation_results = []

        for p_ac in transcripts:
            errors = []

            # Validate aa0 exists at pos0 on given
            invalid_aa0_seq_msg = self.validate_reference_sequence(
                p_ac, classification.pos0, classification.pos0, classification.aa0
            )
            if invalid_aa0_seq_msg:
                errors.append(invalid_aa0_seq_msg)

            # Validate aa1 exists at pos1
            if classification.aa1:
                # TODO: There shouldn't be a case where aa1 exists and pos1 doesn't
                # but we should double check this
                invalid_aa1_seq_msg = self.validate_reference_sequence(
                    p_ac, classification.pos1, classification.pos1, classification.aa1
                )

                if invalid_aa1_seq_msg:
                    errors.append(invalid_aa1_seq_msg)

            validation_results.append(
                ValidationResult(
                    accession=p_ac,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors,
                    gene_tokens=gene_tokens
                )
            )

        return validation_results

    def variation_name(self) -> str:
        """Return the variation name."""
        return "protein insertion"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is protein insertion."""
        return classification_type == ClassificationType.PROTEIN_INSERTION

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
            transcripts = self.get_protein_transcripts(gene_tokens, errors)
        return transcripts

    def get_gene_tokens(self, classification: Classification) -> List:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_protein_gene_symbol_tokens(classification)
