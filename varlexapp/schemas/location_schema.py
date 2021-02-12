"""A module for the Location schema."""
from marshmallow import Schema, fields
from .simple_interval_schema import SimpleIntervalSchema


class LocationSchema(Schema):
    """The location schema class."""

    interval = fields.Nested(SimpleIntervalSchema)
    sequence_id = fields.Str()
    type = fields.Str()
