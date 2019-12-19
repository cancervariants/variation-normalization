from marshmallow import Schema, fields
from .basic_token_schema import BasicTokenSchema
from .gene_match_token_schema import GeneMatchTokenSchema

class GenePairMatchTokenSchema(BasicTokenSchema):
    leftGeneToken = fields.Nested(GeneMatchTokenSchema, attribute='left_gene_token', description='The first of the recognized genes in the pair')
    rightGeneToken = fields.Nested(GeneMatchTokenSchema, attribute='right_gene_token', description='The second of the recognized genes in the pair')
