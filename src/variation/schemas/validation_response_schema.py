"""Module for Validation Response Schema."""

from pydantic import BaseModel, StrictBool, StrictInt, StrictStr

from variation.schemas.classification_response_schema import Classification


class ValidationResult(BaseModel):
    """Validation Results for a given input"""

    accession: StrictStr | None = None
    cds_start: StrictInt | None = None  # This is only for cDNA
    classification: Classification
    is_valid: StrictBool
    errors: list[StrictStr] = []


class ValidationSummary(BaseModel):
    """Give Valid and Invalid Results for a given input."""

    valid_results: list[ValidationResult] = []
    invalid_results: list[ValidationResult] = []
    warnings: list[StrictStr] = []
