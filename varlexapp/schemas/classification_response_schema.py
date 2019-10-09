from marshmallow import Schema, fields

from .classification_schema import ClassificationSchema

class ClassificationResponseSchema(Schema):
    searchTerm = fields.Str(attribute='search_term')
    classifications = fields.Nested(ClassificationSchema, many=True)
