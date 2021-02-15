"""Module for Classification schema."""
from pydantic import BaseModel
from typing import List, Union
from enum import IntEnum
from varlexapp.schemas.token_response_schema import Token, \
    GeneMatchToken, GenePairMatchToken, ProteinSubstitutionToken


class ClassificationType(IntEnum):
    """Define constraints for Classification Types."""

    FUSION = 1
    PROTEIN_SUBSTITUTION = 2
    PROTEIN_FRAMESHIFT = 3
    PROTEIN_ALTERNATE = 4
    PROTEIN_DELINS = 5
    PROTEIN_TERMINATION = 6
    PROTEIN_DUPLICATION = 7
    ONCOGENIC = 8
    EXPRESSION = 9
    COMPLEX = 10


class ConfidenceRating(IntEnum):
    """Define constraints for confidence ratings."""

    INTERSECTION = 1
    SUPERSET = 2
    UNORDERED = 3
    EXACT = 4


class Classification(BaseModel):
    """The classification schema class."""

    classification_type: ClassificationType
    matching_tokens: List[str]
    non_matching_tokens: List[str]
    all_tokens: List[Union[Token, GeneMatchToken, GenePairMatchToken,
                     ProteinSubstitutionToken]]
    confidence: ConfidenceRating

    class Config:
        """Configure model."""

        orm_mode = True


class ClassificationResponseSchema(BaseModel):
    """Define classification response schema."""

    search_term: str
    classifications: List[Classification]

    class Config:
        """Configure model."""

        orm_mode = True
