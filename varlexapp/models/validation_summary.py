from typing import List

from . import ValidationResult

class ValidationSummary:
    def __init__(self, valid_results: List[ValidationResult] , invalid_results: List[ValidationResult]) -> None:
        self.valid_results = valid_results
        self.invalid_results = invalid_results
