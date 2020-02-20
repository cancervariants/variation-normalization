from marshmallow import Schema, fields
from ..models import ConfidenceRating, ClassificationType


from .token_schema import TokenSchema

from .apispec_enum_field import ApispecEnumField

class ClassificationSchema(Schema):
    classificationType = ApispecEnumField(ClassificationType, attribute='classification_type')
    allTokens = fields.Nested(TokenSchema, many=True, attribute='all_tokens')
    matchingTokens = fields.List(fields.Str(), attribute='matching_tokens')
    nonMatchingTokens = fields.List(fields.Str(), attribute='non_matching_tokens')
    confidence = ApispecEnumField(ConfidenceRating)

