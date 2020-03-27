from marshmallow import Schema, fields

from .validation_summary_schema import ValidationSummarySchema

class ValidationResponseSchema(Schema):
    searchTerm = fields.Str(attribute='search_term')
    validationSummary = fields.Nested(ValidationSummarySchema, attribute='validation_summary')
