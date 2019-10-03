from marshmallow import Schema, fields

class ClassificationSchema(Schema):
    classificationType = fields.Str(attribute='classification_type')
    matchingTokens = fields.List(fields.Str(), attribute='matching_tokens')
    nonMatchingTokens = fields.List(fields.Str(), attribute='non_matching_tokens')
    confidence = fields.Str()

