"""Module for classification schema."""
from pydantic import BaseModel
from typing import List, Union
from varlexapp.models import ClassificationType
from varlexapp.schemas.token_response_schema import Token, \
    GeneMatchToken, GenePairMatchToken, ProteinSubstitutionToken


class ClassificationSchema(BaseModel):
    """The classification schema class."""

    classificationType: ClassificationType
    allTokens: List[Union[Token, GeneMatchToken, GenePairMatchToken,
                    ProteinSubstitutionToken]]

    class Config:
        """Configure model."""

        orm_mode = True
