"""Module for modeling VRS objects."""
from pydantic import BaseModel
from pydantic.fields import Field
from pydantic.types import StrictStr
from typing import Dict, Any, Type, Optional
from gene.schemas import SimpleInterval, SequenceLocation  # noqa: F401


class SequenceState(BaseModel):
    """Captures a Sequence as a State."""

    sequence: StrictStr
    type = 'SequenceState'

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['SequenceState']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "sequence": "T",
                "type": "SequenceState"
            }


class Allele(BaseModel):
    """A Sequence or Sequence change with respect to a reference sequence, without regard to genes or other features."""  # noqa: E501

    id: Optional[StrictStr] = Field(alias='_id')
    location: SequenceLocation
    state: SequenceState
    type = 'Allele'

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['Allele']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "location": {
                    "interval": {
                        "end": 44908822,
                        "start": 44908821,
                        "type": "SimpleInterval"
                    },
                    "sequence_id": "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                    "type": "SequenceLocation"
                },
                "state": {
                    "sequence": "T",
                    "type": "SequenceState"
                },
                "type": "Allele"
            }


class Text(BaseModel):
    """A free-text description of variation that is intended for
    interpretation by humans.
    """

    id: Optional[str] = Field(alias='_id')
    type = "Text"
    definition: str
