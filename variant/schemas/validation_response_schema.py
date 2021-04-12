"""Module for Validation Response Schema."""
from pydantic import BaseModel
from pydantic.types import StrictBool
from typing import List, Optional
from enum import IntEnum
from variant.schemas.classification_response_schema import Classification
from variant.schemas.token_response_schema import GeneMatchToken


class LookupType(IntEnum):
    """IntEnum for Lookup Type."""

    GENE_SYMBOL = 1


class ValidationResult(BaseModel):
    """Validation Results for a given variant."""

    classification: Classification
    is_valid: StrictBool
    confidence_score: float
    allele: Optional[dict] = None
    human_description: Optional[str]
    concise_description: str
    errors: List[str]
    gene_tokens: Optional[List[GeneMatchToken]]
    mane_transcript: Optional[str]


class ValidationSummary(BaseModel):
    """Give Valid and Invalid Results for a given variant."""

    valid_results: List[ValidationResult]
    invalid_results: List[ValidationResult]


class ValidationResponseSchema(BaseModel):
    """Validation Response for a given variant."""

    search_term: str
    validation_summary: ValidationSummary
