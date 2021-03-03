"""Module for Validation Response Schema."""
from pydantic import BaseModel
from pydantic.types import StrictBool
from typing import List, Optional, Dict, Any, Type
from enum import IntEnum
from variant.schemas.classification_response_schema import Classification
from variant.schemas.token_response_schema import GeneMatchToken


class LookupType(IntEnum):
    """IntEnum for Lookup Type."""

    GENE_SYMBOL = 1


class ValidationResult(BaseModel):
    """Validation Results for a given variant."""

    classification: Classification
    is_valid: StrictBool
    confidence_score: float
    allele: Optional[dict] = None
    human_description: Optional[str]
    concise_description: str
    errors: List[str]
    gene_tokens: Optional[List[GeneMatchToken]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['ValidationResult']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "classification": {
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
                },
                "is_valid": True,
                "confidence_score": 1,
                "allele": {
                    "_id": "ga4gh:VA.u6sKlz0mMQvARmrlnt0Aksz6EbSkmL8z",
                    "location": {
                        "interval": {
                            "end": 600,
                            "start": 599,
                            "type": "SimpleInterval"
                        },
                        "sequence_id":
                        "ga4gh:SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",
                        "type": "SequenceLocation"
                    },
                    "state": {
                        "sequence": "E",
                        "type": "SequenceState"
                    },
                    "type": "Allele"
                },
                "human_description": "ENSP00000419060.2 V600E",
                "concise_description": "ENSP00000419060.2 V600E",
                "errors": []
            }


class ValidationSummary(BaseModel):
    """Give Valid and Invalid Results for a given variant."""

    valid_results: List[ValidationResult]
    invalid_results: List[ValidationResult]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['ValidationSummary']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "valid_results": [
                    {
                        "classification": {
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
                        },
                        "is_valid": True,
                        "confidence_score": 1,
                        "allele": {
                            "_id": "ga4gh:VA.u6sKlz0mMQvARmrlnt0Aksz6EbSkmL8z",
                            "location": {
                                "interval": {
                                    "end": 600,
                                    "start": 599,
                                    "type": "SimpleInterval"
                                },
                                "sequence_id": "ga4gh:SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",  # noqa: E501
                                "type": "SequenceLocation"
                            },
                            "state": {
                                "sequence": "E",
                                "type": "SequenceState"
                            },
                            "type": "Allele"
                        },
                        "human_description": "ENSP00000419060.2 V600E",
                        "concise_description": "ENSP00000419060.2 V600E",
                        "errors": []
                    }
                ],
                "invalid_results": [
                    {
                        "classification": {
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
                        },
                        "is_valid": False,
                        "confidence_score": 1,
                        "allele": {
                            "_id": "ga4gh:VA.yj1Z7SReNgZx4m6nJ1rfJv3hXglXFwFc",
                            "location": {
                                "interval": {
                                    "end": 600,
                                    "start": 599,
                                    "type": "SimpleInterval"
                                },
                                "sequence_id": "ga4gh:SQ.WaAJ_cXXn9YpMNfhcq9lnzIvaB9ALawo",  # noqa: E501
                                "type": "SequenceLocation"
                            },
                            "state": {
                                "sequence": "E",
                                "type": "SequenceState"
                            },
                            "type": "Allele"
                        },
                        "human_description": "ENSP00000496776.1 V600E",
                        "concise_description": "ENSP00000496776.1 V600E",
                        "errors": [
                            "Needed to find V at position 600 on ENSP00000496776.1 but found T"  # noqa: E501
                        ]
                    }
                ]
            }


class ValidationResponseSchema(BaseModel):
    """Validation Response for a given variant."""

    search_term: str
    validation_summary: ValidationSummary

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['ValidationResponseSchema']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "search_term": "BRAF V600E",
                "valid_results": [
                    {
                        "classification": {
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
                        },
                        "is_valid": True,
                        "confidence_score": 1,
                        "allele": {
                            "_id": "ga4gh:VA.u6sKlz0mMQvARmrlnt0Aksz6EbSkmL8z",
                            "location": {
                                "interval": {
                                    "end": 600,
                                    "start": 599,
                                    "type": "SimpleInterval"
                                },
                                "sequence_id": "ga4gh:SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",  # noqa: E501
                                "type": "SequenceLocation"
                            },
                            "state": {
                                "sequence": "E",
                                "type": "SequenceState"
                            },
                            "type": "Allele"
                        },
                        "human_description": "ENSP00000419060.2 V600E",
                        "concise_description": "ENSP00000419060.2 V600E",
                        "errors": []
                    }
                ],
                "invalid_results": [
                    {
                        "classification": {
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
                        },
                        "is_valid": False,
                        "confidence_score": 1,
                        "allele": {
                            "_id": "ga4gh:VA.yj1Z7SReNgZx4m6nJ1rfJv3hXglXFwFc",
                            "location": {
                                "interval": {
                                    "end": 600,
                                    "start": 599,
                                    "type": "SimpleInterval"
                                },
                                "sequence_id": "ga4gh:SQ.WaAJ_cXXn9YpMNfhcq9lnzIvaB9ALawo",  # noqa: E501
                                "type": "SequenceLocation"
                            },
                            "state": {
                                "sequence": "E",
                                "type": "SequenceState"
                            },
                            "type": "Allele"
                        },
                        "human_description": "ENSP00000496776.1 V600E",
                        "concise_description": "ENSP00000496776.1 V600E",
                        "errors": [
                            "Needed to find V at position 600 on ENSP00000496776.1 but found T"  # noqa: E501
                        ]
                    }
                ]
            }
