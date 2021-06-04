"""Module for VRS location schemas."""
from pydantic import BaseModel
from pydantic.fields import Field
from typing import Dict, Any, Type


class SequenceState(BaseModel):
    """Captures a Sequence as a State."""

    sequence: str
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


class SimpleInterval(BaseModel):
    """A SequenceInterval with a single start and end coordinate."""

    start: int
    end: int
    type = 'SimpleInterval'

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['SimpleInterval']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "end": 44908822,
                "start": 44908821,
                "type": "SimpleInterval"
            }


class SequenceLocation(BaseModel):
    """A specified subsequence within another sequence that is used as a reference sequence."""  # noqa: E501

    interval: SimpleInterval
    sequence_id: str
    type = 'SequenceLocation'

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['SequenceLocation']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "interval": {
                    "end": 44908822,
                    "start": 44908821,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                "type": "SequenceLocation"
            }


class Allele(BaseModel):
    """A Sequence or Sequence change with respect to a reference sequence, without regard to genes or other features."""  # noqa: E501

    id: str = Field(..., alias='_id')
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
