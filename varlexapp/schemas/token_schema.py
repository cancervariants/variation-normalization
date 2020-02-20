from marshmallow_oneofschema import OneOfSchema

from .basic_token_schema import BasicTokenSchema
from .gene_match_token_schema import GeneMatchTokenSchema
from .gene_pair_match_token_schema import GenePairMatchTokenSchema
from .protein_substitution_token_schema import ProteinSubstitutionTokenSchema


class TokenSchema(OneOfSchema):
    type_schemas = {
            'BasicToken': BasicTokenSchema,
            'GeneMatchToken': GeneMatchTokenSchema,
            'GenePairMatchToken': GenePairMatchTokenSchema,
            'ProteinSubstitutionToken': ProteinSubstitutionTokenSchema
        }
