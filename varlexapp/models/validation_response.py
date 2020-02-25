from . import ValidationSummary

class ValidationResponse:
    def __init__(self, search_term: str, validation_summary: ValidationSummary) -> None:
        self.search_term = search_term
        self.validation_summary = validation_summary
