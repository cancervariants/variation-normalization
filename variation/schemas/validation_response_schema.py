"""Module for Validation Response Schema."""
from typing import List, Optional
from enum import IntEnum

from pydantic import BaseModel
from pydantic.types import StrictBool

from variation.schemas.classification_response_schema import Classification
from variation.schemas.token_response_schema import GeneMatchToken, Token


class LookupType(IntEnum):
    """IntEnum for Lookup Type."""

    GENE_SYMBOL = 1


class ValidationResult(BaseModel):
    """Validation Results for a given variation."""

    classification: Classification
    classification_token: Optional[Token]
    is_valid: StrictBool
    confidence_score: float
    variation: Optional[dict] = None
    errors: List[str]
    gene_tokens: Optional[List[GeneMatchToken]]
    is_mane_transcript: Optional[StrictBool]
    identifier: Optional[str]


class ValidationSummary(BaseModel):
    """Give Valid and Invalid Results for a given variation."""

    valid_results: List[ValidationResult]
    invalid_results: List[ValidationResult]
    warnings: Optional[List[str]]


class ValidationResponseSchema(BaseModel):
    """Validation Response for a given variation."""

    search_term: str
    validation_summary: ValidationSummary
