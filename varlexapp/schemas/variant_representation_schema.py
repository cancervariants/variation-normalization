"""Module for Variant Representation schema."""
from marshmallow import Schema, fields
from .location_schema import LocationSchema
from .sequence_state_schema import SequenceStateSchema


class VariantRepresentationSchema(Schema):
    """The Variant Representation Schema class."""

    _id = fields.Str()
    location = fields.Nested(LocationSchema)
    state = fields.Nested(SequenceStateSchema)
    type = fields.Str()
