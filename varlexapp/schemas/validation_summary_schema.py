from marshmallow import Schema, fields

from .validation_result_schema import ValidationResultSchema

class ValidationSummarySchema(Schema):
    validResults = fields.Nested(ValidationResultSchema, many=True, attribute='valid_results')
    invalidResults = fields.Nested(ValidationResultSchema, many=True, attribute='invalid_results')



