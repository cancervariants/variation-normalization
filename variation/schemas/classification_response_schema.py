"""Module for Classification schema."""
from typing import List
from enum import IntEnum, Enum

from pydantic import BaseModel

from variation.schemas.token_response_schema import Token


class ClassificationType(str, Enum):
    """Enums for Classification Types."""

    FUSION = "fusion"
    PROTEIN_SUBSTITUTION = "protein substitution"
    POLYPEPTIDE_TRUNCATION = "polypeptide truncation"
    SILENT_MUTATION = "silent mutation"
    PROTEIN_FRAMESHIFT = "protein frameshift"
    PROTEIN_ALTERNATE = "protein alternate"
    PROTEIN_DELINS = "protein delins"
    PROTEIN_TERMINATION = "protein termination"
    PROTEIN_DUPLICATION = "protein duplication"
    ONCOGENIC = "oncogenic"
    EXPRESSION = "expression"
    COMPLEX = "complext"
    CODING_DNA_SUBSTITUTION = "coding dna substitution"
    GENOMIC_SUBSTITUTION = "genomic substitution"
    CODING_DNA_SILENT_MUTATION = "coding dna silent mutation"
    GENOMIC_SILENT_MUTATION = "genomic silent mutation"
    CODING_DNA_DELINS = "coding dna delins"
    GENOMIC_DELINS = "genomic delins"
    PROTEIN_DELETION = "protein deletion"
    CODING_DNA_DELETION = "coding dna deletion"
    GENOMIC_DELETION = "genomic deletion"
    PROTEIN_INSERTION = "protein insertion"
    CODING_DNA_INSERTION = "coding dna insertion"
    GENOMIC_INSERTION = "genomic insertion"
    GENOMIC_UNCERTAIN_DELETION = "genomic uncertain deletion"
    GENOMIC_DUPLICATION = "genomic duplication"
    GENOMIC_DELETION_RANGE = "genomic deletion range"


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
