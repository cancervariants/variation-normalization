from marshmallow import Schema, fields

from .simple_interval_schema import SimpleIntervalSchema

class LocationSchema(Schema):
    interval = fields.Nested(SimpleIntervalSchema)
    sequence_id = fields.Str()
    type = fields.Str()

