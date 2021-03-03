"""Module for Classification schema."""
from pydantic import BaseModel
from typing import List, Union, Dict, Any, Type
from enum import IntEnum
from variant.schemas.token_response_schema import Token, \
    GeneMatchToken, GenePairMatchToken, AminoAcidSubstitutionToken, \
    PolypeptideTruncationToken, SilentMutationToken


class ClassificationType(IntEnum):
    """Enums for Classification Types."""

    FUSION = 1
    PROTEIN_SUBSTITUTION = 2
    AMINO_ACID_SUBSTITUTION = 3
    POLYPEPTIDE_TRUNCATION = 4
    SILENT_MUTATION = 5
    PROTEIN_FRAMESHIFT = 6
    PROTEIN_ALTERNATE = 7
    PROTEIN_DELINS = 8
    PROTEIN_TERMINATION = 9
    PROTEIN_DUPLICATION = 10
    ONCOGENIC = 11
    EXPRESSION = 12
    COMPLEX = 13


class ConfidenceRating(IntEnum):
    """Enums for classification confidence ratings."""

    INTERSECTION = 1
    SUPERSET = 2
    UNORDERED = 3
    EXACT = 4


class Classification(BaseModel):
    """Classification for a list of tokens."""

    classification_type: ClassificationType
    matching_tokens: List[str]
    non_matching_tokens: List[str]
    all_tokens: List[Union[GeneMatchToken, GenePairMatchToken,
                           AminoAcidSubstitutionToken,
                           PolypeptideTruncationToken,
                           SilentMutationToken, Token]]
    confidence: ConfidenceRating

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['Classification']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "classification_type": 2,
                "matching_tokens": [
                    "GeneSymbol",
                    "ProteinSubstitution"
                ],
                "non_matching_tokens": [],
                "all_tokens": [
                    {
                        "token": "BRAF",
                        "token_type": "GeneSymbol",
                        "match_type": 2,
                        "input_string": "BRAF",
                        "object_type": "Token",
                        "matched_value": "BRAF"
                    },
                    {
                        "token": "V600E",
                        "token_type": "ProteinSubstitution",
                        "match_type": 5,
                        "input_string": "V600E",
                        "object_type": "Token",
                        "ref_protein": "V",
                        "alt_protein": "E",
                        "position": 600
                    }
                ],
                "confidence": 4
            }


class ClassificationResponseSchema(BaseModel):
    """Classification response for a given query."""

    search_term: str
    classifications: List[Classification]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['ClassificationResponseSchema']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "search_term": "BRAF V600E",
                "classifications": [
                    {
                        "classification_type": 2,
                        "matching_tokens": [
                            "GeneSymbol",
                            "ProteinSubstitution"
                        ],
                        "non_matching_tokens": [],
                        "all_tokens": [
                            {
                                "token": "BRAF",
                                "token_type": "GeneSymbol",
                                "match_type": 2,
                                "input_string": "BRAF",
                                "object_type": "Token",
                                "matched_value": "BRAF"
                            },
                            {
                                "token": "V600E",
                                "token_type": "ProteinSubstitution",
                                "match_type": 5,
                                "input_string": "V600E",
                                "object_type": "Token"
                            }
                        ],
                        "confidence": 4
                    }
                ]
            }
