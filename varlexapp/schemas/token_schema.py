from marshmallow import Schema, fields

class TokenSchema(Schema):
    token = fields.Str()
    tokenType = fields.Str(attribute='token_type')
