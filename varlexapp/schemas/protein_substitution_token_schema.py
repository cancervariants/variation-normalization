from marshmallow import Schema, fields
from .basic_token_schema import BasicTokenSchema

class ProteinSubstitutionTokenSchema(BasicTokenSchema):
    refProtein = fields.Str(attribute='ref_protein', description='Expected protein at pos')
    altProtein = fields.Str(attribute='alt_protein', description='Protein that replaces ref_protein')
    pos = fields.Str(description='Position ref_protein is expected')
