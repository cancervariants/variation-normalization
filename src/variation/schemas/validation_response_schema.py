"""Module for Validation Response Schema."""
from typing import List, Optional

from pydantic import BaseModel, StrictBool, StrictInt, StrictStr

from variation.schemas.classification_response_schema import Classification


class ValidationResult(BaseModel):
    """Validation Results for a given input"""

    accession: Optional[StrictStr] = None
    cds_start: Optional[StrictInt] = None  # This is only for cDNA
    classification: Classification
    is_valid: StrictBool
    errors: List[StrictStr] = []


class ValidationSummary(BaseModel):
    """Give Valid and Invalid Results for a given input."""

    valid_results: List[ValidationResult] = []
    invalid_results: List[ValidationResult] = []
    warnings: List[StrictStr] = []
