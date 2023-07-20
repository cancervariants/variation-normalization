"""The module for Protein Stop Gain Validation."""
from typing import Optional, List

from variation.schemas.classification_response_schema import (
    ClassificationType, Classification, Nomenclature, ProteinStopGainClassification
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.utils import get_aa1_codes
from variation.validators.validator import Validator


class ProteinStopGain(Validator):
    """The Protein Stop Gain Validator class."""

    async def get_valid_invalid_results(
        self, classification: ProteinStopGainClassification, transcripts: List[str]
    ) -> List[ValidationResult]:
        errors = []

        # Only HGVS Expressions are validated
        # Free text is validated during tokenization
        # Don't need to validate alt, since we know it's '*'
        if classification.nomenclature == Nomenclature.HGVS:
            aa1_ref = get_aa1_codes(classification.ref)
            if aa1_ref:
                classification.ref = aa1_ref
            else:
                errors.append(f"`ref` not valid amino acid(s): {classification.ref}")

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
        return "protein stop gain"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is protein stop gain."""
        return classification_type == ClassificationType.PROTEIN_STOP_GAIN

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
