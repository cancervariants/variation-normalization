"""Module for Validation Response Schema."""
from pydantic import BaseModel
from pydantic.types import StrictBool
from typing import List, Optional
from enum import IntEnum
from varlexapp.schemas.classification_response_schema import Classification


class LookupType(IntEnum):
    """Define constraints for Lookup Type."""

    GENE_SYMBOL = 1


class ValidationResult(BaseModel):
    """Define model for validation result."""

    classification: Classification
    is_valid: StrictBool
    confidence_score: float
    allele: Optional[dict] = None
    human_description: Optional[str]
    concise_description: str
    errors: List[str]

    class Config:
        """Configure model."""

        orm_mode = True


class ValidationSummary(BaseModel):
    """Define model for validation summary."""

    valid_results: List[ValidationResult]
    invalid_results: List[ValidationResult]

    class Config:
        """Configure model."""

        orm_mode = True


class ValidationResponseSchema(BaseModel):
    """Define model for validation response."""

    search_term: str
    validation_summary: ValidationSummary

    class Config:
        """Configure model."""

        orm_mode = True
