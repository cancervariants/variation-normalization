"""Module for VRS location schemas."""
from pydantic import BaseModel
from pydantic.fields import Field


class SequenceState(BaseModel):
    """VRS Sequence State constraints."""

    sequence: str
    type = 'SequenceState'

    class Config:
        """Configure model."""

        orm_mode = True


class SimpleInterval(BaseModel):
    """VRS Simple Interval constraints."""

    start: int
    end: int
    type = 'SimpleInterval'

    class Config:
        """Configure model."""

        orm_mode = True


class SequenceLocation(BaseModel):
    """VRS sequence location constraints."""

    interval: SimpleInterval
    sequence_id: str
    type = 'SequenceLocation'

    class Config:
        """Configure model."""

        orm_mode = True


class Allele(BaseModel):
    """VRS Allele constraints."""

    id: str = Field(..., alias='_id')
    location: SequenceLocation
    state: SequenceState
    type = 'Allele'

    class Config:
        """Configure model."""

        orm_mode = True
