from marshmallow import Schema, fields
from ..models import ConfidenceRating, ClassificationType

from .apispec_enum_field import ApispecEnumField

class ClassificationSchema(Schema):
    classificationType = ApispecEnumField(ClassificationType, attribute='classification_type')
    matchingTokens = fields.List(fields.Str(), attribute='matching_tokens')
    nonMatchingTokens = fields.List(fields.Str(), attribute='non_matching_tokens')
    confidence = ApispecEnumField(ConfidenceRating)

