"""Module for Classification schema."""
from typing import List
from enum import IntEnum

from pydantic import BaseModel

from variation.schemas.token_response_schema import Token


class ClassificationType(IntEnum):
    """Enums for Classification Types."""

    FUSION = 1
    PROTEIN_SUBSTITUTION = 2
    POLYPEPTIDE_TRUNCATION = 3
    SILENT_MUTATION = 4
    PROTEIN_FRAMESHIFT = 5
    PROTEIN_ALTERNATE = 6
    PROTEIN_DELINS = 7
    PROTEIN_TERMINATION = 8
    PROTEIN_DUPLICATION = 9
    ONCOGENIC = 10
    EXPRESSION = 11
    COMPLEX = 12
    CODING_DNA_SUBSTITUTION = 13
    GENOMIC_SUBSTITUTION = 14
    CODING_DNA_SILENT_MUTATION = 15
    GENOMIC_SILENT_MUTATION = 16
    CODING_DNA_DELINS = 17
    GENOMIC_DELINS = 18
    PROTEIN_DELETION = 19
    CODING_DNA_DELETION = 20
    GENOMIC_DELETION = 21
    PROTEIN_INSERTION = 22
    CODING_DNA_INSERTION = 23
    GENOMIC_INSERTION = 24
    GENOMIC_UNCERTAIN_DELETION = 25
    GENOMIC_DUPLICATION = 26
    GENOMIC_DELETION_RANGE = 27


class ConfidenceRating(IntEnum):
    """Enums for classification confidence ratings."""

    INTERSECTION = 1
    SUPERSET = 2
    UNORDERED = 3
    EXACT = 4


class Classification(BaseModel):
    """Classification for a list of tokens."""

    classification_type: ClassificationType
    matching_tokens: List[str]
    non_matching_tokens: List[str]
    all_tokens: List[Token]
    confidence: ConfidenceRating


class ClassificationResponseSchema(BaseModel):
    """Classification response for a given query."""

    search_term: str
    classifications: List[Classification]
