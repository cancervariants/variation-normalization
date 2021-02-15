"""Module for validate result."""
from marshmallow import Schema, fields
from varlexapp.schemas.classification_response_schema import Classification


class ValidationResultSchema(Schema):
    """The validation result schema class."""

    classification = fields.Nested(Classification)
    isValid = fields.Boolean(attribute='is_valid')
    confidenceScore = fields.Float(attribute='confidence_score')
    humanDescription = fields.Str(attribute='human_description')
    conciseDescription = fields.Str(attribute='concise_description')
    errors = fields.List(fields.Str())
