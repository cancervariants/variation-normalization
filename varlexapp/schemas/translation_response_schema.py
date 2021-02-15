"""Module for Translation Response Schema."""
from pydantic import BaseModel
from typing import List
from varlexapp.schemas.ga4gh_vrs import Allele


class TranslationResponseSchema(BaseModel):
    """Define model for translation response."""

    search_term: str
    variants: List[Allele]

    class Config:
        """Configure model."""

        orm_mode = True
