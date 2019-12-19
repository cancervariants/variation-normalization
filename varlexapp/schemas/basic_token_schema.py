from marshmallow import Schema, fields
from ..models import TokenMatchType
from .apispec_enum_field import ApispecEnumField

class BasicTokenSchema(Schema):
    token = fields.Str(description='The token identified for this string')
    tokenType = fields.Str(description='The type of token recognized', attribute='token_type')
    matchType = ApispecEnumField(TokenMatchType, description='The type of comparison used to determine the match', attribute='match_type')
    inputString = fields.Str(description='The original string', attribute='input_string')
    objectType = fields.Str(description='What type of token is this', attribute='object_type')
