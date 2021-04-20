"""Module for Token Schema."""
from pydantic import BaseModel
from typing import List, Union, Dict, Any, Type, Optional
from enum import IntEnum, Enum


class TokenMatchType(IntEnum):
    """Gene Token match types."""

    ID = 1
    SYMBOL = 2
    ALIAS = 3
    PREVIOUS = 4
    UNSPECIFIED = 5


class Token(BaseModel):
    """A string from a given query."""

    token: str
    token_type: str
    match_type: TokenMatchType
    input_string: str
    object_type = 'Token'

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['Token']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "token": "",
                "token_type": "Unknown",
                "match_type": 5,
                "input_string": "foo",
                "object_type": "Token"
            }


class GeneMatchToken(Token):
    """Define model for gene symbol token."""

    matched_value: str
    token_type = 'GeneSymbol'

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['GeneMatchToken']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "token": "BRAF",
                "token_type": "GeneSymbol",
                "match_type": 2,
                "input_string": "BRAF",
                "object_type": "Token",
                "matched_value": "BRAF"
            }


class GenePairMatchToken(Token):
    """Define model for gene pair token."""

    left_gene_token: GeneMatchToken
    right_gene_token: GeneMatchToken
    token_type = 'GenePair'

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['GenePairMatchToken']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "token": "BRAF-ABL1",
                "token_type": "GenePair",
                "match_type": 5,
                "input_string": "braf-abl1",
                "object_type": "Token",
                "left_gene_token": {
                    "token": "BRAF",
                    "token_type": "GeneSymbol",
                    "match_type": 2,
                    "input_string": "BRAF",
                    "object_type": "Token",
                    "matched_value": "BRAF"
                },
                "right_gene_token": {
                    "token": "ABL1",
                    "token_type": "GeneSymbol",
                    "match_type": 2,
                    "input_string": "ABL1",
                    "object_type": "Token",
                    "matched_value": "ABL1"
                }
            }


class PolypeptideSequenceVariant(Token):
    """Polypeptide Sequence Variant Token Class."""

    ref_protein: str
    alt_protein: str
    position: int
    token_type: str


class PolypeptideTruncationToken(PolypeptideSequenceVariant):
    """A sequence variant of the CD that causes a truncation of the
    resulting polypeptide. (nonsense)
    """

    alt_protein = 'Ter'
    token_type = 'PolypeptideTruncation'

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['PolypeptideTruncationToken']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "token": "Tyr365Ter",
                "token_type": "PolypeptideTruncation",
                "match_type": 5,
                "input_string": "Tyr365Ter",
                "object_type": "Token",
                "ref_protein": "Tyr",
                "alt_protein": "Ter",
                "position": 365
            }


class AminoAcidSubstitutionToken(PolypeptideSequenceVariant):
    """A sequence variant of a codon resulting in the substitution of one
    amino acid for another in the resulting polypeptide. (missense)
    """

    token_type = 'AminoAcidSubstitution'

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['AminoAcidSubstitutionToken']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "token": "V600E",
                "token_type": "AminoAcidSubstitution",
                "match_type": 5,
                "input_string": "V600E",
                "object_type": "Token",
                "ref_protein": "V",
                "alt_protein": "E",
                "position": 600
            }


class SilentMutationToken(PolypeptideSequenceVariant):
    """A sequence variant that does not affect protein functions."""

    alt_protein = '='
    token_type = 'SilentMutation'

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['SilentMutationToken']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "token": "p.Cys188=",
                "token_type": "SilentMutation",
                "match_type": 5,
                "input_string": "p.Cys188=",
                "object_type": "Token",
                "ref_protein": "Cys",
                "alt_protein": "=",
                "position": 188
            }


class TokenResponseSchema(BaseModel):
    """Define model for token response."""

    search_term: str
    tokens: List[Union[GeneMatchToken, GenePairMatchToken,
                       AminoAcidSubstitutionToken, PolypeptideTruncationToken,
                       SilentMutationToken, Token]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['TokenResponseSchema']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "search_term": "BRAF V600E",
                "tokens": [
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
                        "token_type": "AminoAcidSubstitution",
                        "match_type": 5,
                        "input_string": "V600E",
                        "object_type": "Token",
                        "ref_protein": "V",
                        "alt_protein": "E",
                        "position": 600
                    }
                ]
            }


class ReferenceSequence(str, Enum):
    """Define constraints for reference sequence."""

    CODING_DNA = 'c'
    LINEAR_GENOMIC = 'g'


class SequenceAlteration(Token):
    """Sequence feature whose extent is the deviation from another sequence.
    (SO:0001059)
    """

    position: int
    ref_nucleotide: Optional[str]
    new_nucleotide: str
    token_type: str
    reference_sequence: ReferenceSequence


class SingleNucleotideVariantToken(SequenceAlteration):
    """Single nucleotide positions in genomic DNA at which different
    sequence alternatives exist.
    """

    token_type = 'SNV'


class CodingDNASubstitutionToken(SingleNucleotideVariantToken):
    """SNV substitution at the coding DNA reference sequence."""

    reference_sequence = ReferenceSequence.CODING_DNA
    token_type = 'CodingDNASubstitution'


class CodingDNASilentMutationToken(SingleNucleotideVariantToken):
    """SNV no change at the coding DNA reference sequence."""

    reference_sequence = ReferenceSequence.CODING_DNA
    new_nucleotide = '='
    token_type = 'CodingDNASilentMutation'


class GenomicSubstitutionToken(SingleNucleotideVariantToken):
    """SNV substitution at the linear genomic reference sequence."""

    reference_sequence = ReferenceSequence.LINEAR_GENOMIC
    token_type = 'GenomicSubstitution'


class GenomicSilentMutationToken(SingleNucleotideVariantToken):
    """SNV no change at the linear genomic reference sequence."""

    reference_sequence = ReferenceSequence.LINEAR_GENOMIC
    new_nucleotide = '='
    token_type = 'GenomicSilentMutation'
