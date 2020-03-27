from marshmallow import Schema, fields

from .variant_representation_schema import VariantRepresentationSchema

class TranslationResponseSchema(Schema):
    searchTerm = fields.Str(attribute='search_term')
    variants = fields.Nested(VariantRepresentationSchema, many=True)
