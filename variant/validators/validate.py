"""Module for Validation."""
from variant.schemas.validation_response_schema import ValidationSummary
from variant.schemas.classification_response_schema import Classification
from ..data_sources import TranscriptMappings, SeqRepoAccess
from .amino_acid_substitution import AminoAcidSubstitution
from .polypeptide_truncation import PolypeptideTruncation
from .silent_mutation import SilentMutation
from typing import List


class Validate:
    """The validation class."""

    def __init__(self, seq_repo_client: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings) -> None:
        """Initialize the validate class."""
        self.validators = [
            AminoAcidSubstitution(seq_repo_client, transcript_mappings),
            PolypeptideTruncation(seq_repo_client, transcript_mappings),
            SilentMutation(seq_repo_client, transcript_mappings)
        ]

    def perform(self, classifications: List[Classification]) \
            -> ValidationSummary:
        """Validate a list of classifications."""
        valid_possibilities = list()
        invalid_possibilities = list()

        for classification in classifications:
            for validator in self.validators:
                if validator.validates_classification_type(
                        classification.classification_type):
                    results = validator.validate(classification)
                    for res in results:
                        if res.is_valid:
                            valid_possibilities.append(res)
                        else:
                            invalid_possibilities.append(res)

        return ValidationSummary(
            valid_results=valid_possibilities,
            invalid_results=invalid_possibilities
        )
