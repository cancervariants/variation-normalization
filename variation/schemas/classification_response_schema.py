"""Module for Classification schema."""
from pydantic import BaseModel
from typing import List
from enum import IntEnum
from variation.schemas.token_response_schema import Token


class ClassificationType(IntEnum):
    """Enums for Classification Types."""

    FUSION = 1
    PROTEIN_SUBSTITUTION = 2
    AMINO_ACID_SUBSTITUTION = 3
    POLYPEPTIDE_TRUNCATION = 4
    SILENT_MUTATION = 5
    PROTEIN_FRAMESHIFT = 6
    PROTEIN_ALTERNATE = 7
    PROTEIN_DELINS = 8
    PROTEIN_TERMINATION = 9
    PROTEIN_DUPLICATION = 10
    ONCOGENIC = 11
    EXPRESSION = 12
    COMPLEX = 13
    CODING_DNA_SUBSTITUTION = 14
    GENOMIC_SUBSTITUTION = 15
    CODING_DNA_SILENT_MUTATION = 16
    GENOMIC_SILENT_MUTATION = 17
    AMINO_ACID_DELINS = 18
    CODING_DNA_DELINS = 19
    GENOMIC_DELINS = 20
    AMINO_ACID_DELETION = 21
    CODING_DNA_DELETION = 22
    GENOMIC_DELETION = 23
    AMINO_ACID_INSERTION = 24
    CODING_DNA_INSERTION = 25
    GENOMIC_INSERTION = 26
    GENOMIC_UNCERTAIN_DELETION = 27


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
