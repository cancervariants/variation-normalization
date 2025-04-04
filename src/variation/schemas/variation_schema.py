"""Define supported variation types"""

from typing import Literal

from pydantic import BaseModel, StrictInt, StrictStr


class Substitution(BaseModel):
    """Define model for substitution variation"""

    pos: StrictInt
    ref: StrictStr
    alt: StrictStr


class StopGain(Substitution):
    """Define model for stop gain variation"""

    alt: Literal["*"] = "*"


class Deletion(BaseModel):
    """Define model for deletion variation"""

    pos0: StrictInt
    pos1: StrictInt | None = None
    deleted_sequence: StrictStr | None = None


class ProteinDeletion(Deletion):
    """Define model for protein deletion"""

    aa0: StrictStr
    aa1: StrictStr | None = None


class Insertion(BaseModel):
    """Define model for insertion variation"""

    pos0: StrictInt
    pos1: StrictInt
    inserted_sequence: StrictStr


class ProteinInsertion(Insertion):
    """Define model for protein insertion variation"""

    aa0: StrictStr
    aa1: StrictStr


class ReferenceAgree(BaseModel):
    """Define model for reference agree variation"""

    pos: StrictInt


class ProteinReferenceAgree(ReferenceAgree):
    """Define model for protein reference agree variation"""

    ref: StrictStr


class DelIns(BaseModel):
    """Define model for delins variation"""

    pos0: StrictInt
    pos1: StrictInt | None = None
    inserted_sequence: StrictStr


class ProteinDelIns(DelIns):
    """Define model for protein delins variation"""

    aa0: StrictStr
    aa1: StrictStr | None = None


class Duplication(BaseModel):
    """Define model for duplication variation"""

    pos0: StrictInt
    pos1: StrictInt | None = None


class DupDelAmbiguous(BaseModel):
    """Define model for duplication/deletion ambiguous variation"""

    pos0: StrictInt | Literal["?"]
    pos1: StrictInt | Literal["?"] | None = None
    pos2: StrictInt | Literal["?"]
    pos3: StrictInt | Literal["?"] | None = None
