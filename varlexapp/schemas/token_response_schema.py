from marshmallow import Schema, fields

from .token_schema import TokenSchema

class TokenResponseSchema(Schema):
    searchTerm = fields.Str(attribute='search_term')
    tokens = fields.Nested(TokenSchema, many=True)
