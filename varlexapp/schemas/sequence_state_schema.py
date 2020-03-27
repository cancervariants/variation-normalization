from marshmallow import Schema, fields


class SequenceStateSchema(Schema):
    sequence = fields.Str()
    type = fields.Str()
