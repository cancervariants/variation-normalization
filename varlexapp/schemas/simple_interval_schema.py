from marshmallow import Schema, fields

class SimpleIntervalSchema(Schema):
    start = fields.Int()
    end = fields.Int()
    type = fields.Str()

