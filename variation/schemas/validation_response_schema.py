"""Module for Validation Response Schema."""
from typing import List, Optional

from pydantic import BaseModel
from pydantic.types import StrictBool, StrictStr, StrictInt

from variation.schemas.classification_response_schema import Classification


class ValidationResult(BaseModel):
    """Validation Results for a given variation."""

    accession: Optional[StrictStr]
    cds_start: Optional[StrictInt]  # This is only for cDNA
    classification: Classification
    is_valid: StrictBool
    errors: List[StrictStr] = []


class ValidationSummary(BaseModel):
    """Give Valid and Invalid Results for a given variation."""

    valid_results: List[ValidationResult]
    invalid_results: List[ValidationResult]
    warnings: List[str]


class ValidationResponseSchema(BaseModel):
    """Validation Response for a given variation."""

    search_term: str
    validation_summary: ValidationSummary
