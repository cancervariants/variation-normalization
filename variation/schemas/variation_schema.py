"""Define supported variation types"""
from typing import Optional, Literal

from pydantic import BaseModel, StrictInt, StrictStr


class Substitution(BaseModel):

    pos: StrictInt
    ref: StrictStr
    alt: StrictStr


class StopGain(Substitution):

    alt: Literal["*"] = "*"


class Deletion(BaseModel):

    pos0: StrictInt
    pos1: Optional[StrictInt]
    deleted_sequence: Optional[StrictStr]


class ProteinDeletion(Deletion):

    aa0: StrictStr
    aa1: Optional[StrictStr]


class Insertion(BaseModel):

    pos0: StrictInt
    pos1: StrictInt
    inserted_sequence: StrictStr


class ProteinInsertion(Insertion):

    aa0: StrictStr
    aa1: StrictStr


class ReferenceAgree(BaseModel):

    pos: StrictInt


class ProteinReferenceAgree(ReferenceAgree):

    ref: StrictStr


class DelIns(BaseModel):

    pos0: StrictInt
    pos1: Optional[StrictInt]
    inserted_sequence: StrictStr


class ProteinDelIns(DelIns):

    aa0: StrictStr
    aa1: Optional[StrictStr]
