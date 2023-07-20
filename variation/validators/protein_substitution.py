"""The module for Protein Substitution Validation."""
from typing import Optional, List

from variation.schemas.classification_response_schema import (
    ClassificationType, Classification, Nomenclature, ProteinSubstitutionClassification
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.utils import get_aa1_codes
from variation.validators.validator import Validator


class ProteinSubstitution(Validator):
    """The Protein Substitution Validator class."""

    async def get_valid_invalid_results(
        self, classification: ProteinSubstitutionClassification,
        transcripts: List[str]
    ) -> List[ValidationResult]:
        errors = []

        # Only HGVS Expressions are validated
        # Free text is validated during tokenization
        if classification.nomenclature == Nomenclature.HGVS:
            aa1_ref = get_aa1_codes(classification.ref)
            if aa1_ref:
                classification.ref = aa1_ref
            else:
                errors.append(f"`ref` not valid amino acid(s): {classification.ref}")

            aa1_alt = get_aa1_codes(classification.alt)
            if aa1_alt:
                classification.alt = aa1_alt
            else:
                errors.append(f"`alt` not valid amino acid(s): {classification.alt}")

            if errors:
                return [ValidationResult(
                    accession=None,
                    classification=classification,
                    is_valid=False,
                    errors=errors
                )]

        validation_results = []

        for t in transcripts:
            errors = []

            valid_ref_seq_msg = self.validate_reference_sequence(
                t, classification.pos, classification.pos, classification.ref
            )
            if valid_ref_seq_msg:
                errors.append(valid_ref_seq_msg)

            validation_results.append(
                ValidationResult(
                    accession=t,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors
                )
            )

        return validation_results

    def variation_name(self) -> str:
        """Return the variation name."""
        return "protein substitution"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is protein substitution."""
        return classification_type == ClassificationType.PROTEIN_SUBSTITUTION

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

