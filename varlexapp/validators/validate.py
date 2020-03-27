from ..models import ValidationResult, ValidationSummary, Classification
from ..data_sources import TranscriptMappings, SeqRepoAccess
from .protein_substitution import ProteinSubstitution

from typing import List

class Validate:
    def __init__(self, seq_repo_client: SeqRepoAccess, transcript_mappings: TranscriptMappings) -> None:
        self.validators = [
                ProteinSubstitution(seq_repo_client, transcript_mappings)
        ]

    def perform(self, classifications: List[Classification]) -> ValidationSummary:
        valid_possibilities = list()
        invalid_possibilities = list()

        for classification in classifications:
            for validator in self.validators:
                if validator.validates_classification_type(classification.classification_type):
                    results = validator.validate(classification)
                    for res in results:
                        if res.is_valid:
                            valid_possibilities.append(res)
                        else:
                            invalid_possibilities.append(res)

        return ValidationSummary(valid_possibilities, invalid_possibilities)

