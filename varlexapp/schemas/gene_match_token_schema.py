from marshmallow import Schema, fields
from .basic_token_schema import BasicTokenSchema

class GeneMatchTokenSchema(BasicTokenSchema):
    matchedValue = fields.Str(description='Test stapling extra metadata in a subtype', attribute='matched_value')
