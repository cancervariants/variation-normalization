"""Module for Token Schema."""
from pydantic import BaseModel
from typing import List, Union, Optional
from enum import IntEnum


class TokenMatchType(IntEnum):
    """Define enum for token match type."""

    ID = 1
    SYMBOL = 2
    ALIAS = 3
    PREVIOUS = 4
    UNSPECIFIED = 5


class Token(BaseModel):
    """Define model for a token."""

    token: str
    token_type: str
    match_type: TokenMatchType
    input_string: str
    object_type = 'Token'

    class Config:
        """Configure model."""

        orm_mode = True


class GeneMatchToken(Token):
    """Define model for gene symbol token."""

    matched_value: str
    token_type = 'GeneSymbol'

    class Config:
        """Configure model."""

        orm_mode = True


class GenePairMatchToken(Token):
    """Define model for gene pair token."""

    left_gene_token: GeneMatchToken
    right_gene_token: GeneMatchToken
    token_type = 'GenePair'

    class Config:
        """Configure model."""

        orm_mode = True


class ProteinSubstitutionToken(Token):
    """Define model for Protein Substitution token."""

    # TODO: These might not be optional
    ref_protein: Optional[str]
    alt_protein: Optional[str]
    position: Optional[int]
    token_type = 'ProteinSubstitution'

    class Config:
        """Configure model."""

        orm_mode = True


class TokenResponseSchema(BaseModel):
    """Define model for token response."""

    search_term: str
    tokens: List[Union[Token, GeneMatchToken, GenePairMatchToken,
                       ProteinSubstitutionToken]]

    class Config:
        """Configure model."""

        orm_mode = True
