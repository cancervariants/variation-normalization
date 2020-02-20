from marshmallow import Schema, fields

from .classification_schema import ClassificationSchema

class ValidationResultSchema(Schema):
    classification = fields.Nested(ClassificationSchema)
    isValid = fields.Boolean(attribute='is_valid')
    confidenceScore = fields.Float(attribute='confidence_score')
    humanDescription = fields.Str(attribute='human_description')
    conciseDescription = fields.Str(attribute='concise_description')
    errors = fields.List(fields.Str())

