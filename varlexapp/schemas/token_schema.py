from marshmallow import Schema, fields
from ..models import TokenMatchType
from .apispec_enum_field import ApispecEnumField

class TokenSchema(Schema):
    token = fields.Str()
    tokenType = fields.Str(attribute='token_type')
    matchType = ApispecEnumField(TokenMatchType, attribute='match_type')
    inputString = fields.Str(attribute='input_string')
