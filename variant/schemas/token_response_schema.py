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


class ReferenceSequence(str, Enum):
    """Define constraints for reference sequence."""

    CODING_DNA = 'c'
    LINEAR_GENOMIC = 'g'
    PROTEIN = 'p'


class PolypeptideSequenceVariant(Token):
    """Polypeptide Sequence Variant Token Class."""

    ref_protein: str
    alt_protein: str
    position: int
    token_type: str
    reference_sequence = ReferenceSequence.PROTEIN
    so_id: str
    molecule_context = 'protein'


class PolypeptideTruncationToken(PolypeptideSequenceVariant):
    """A sequence variant of the CD that causes a truncation of the
    resulting polypeptide. (nonsense)
    """

    alt_protein = 'Ter'
    token_type = 'PolypeptideTruncation'
    so_id = 'SO:0001617'

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
    so_id = 'SO:0001606'

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
    so_id = 'SO:0001017'

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


class SingleNucleotideVariant(Token):
    """Single nucleotide positions in genomic DNA at which different
    sequence alternatives exist.
    """

    position: int
    ref_nucleotide: Optional[str]
    new_nucleotide: str
    token_type: str
    reference_sequence: ReferenceSequence
    so_id: str
    molecule_context: str


class CodingDNASubstitutionToken(SingleNucleotideVariant):
    """SNV substitution at the coding DNA reference sequence."""

    reference_sequence = ReferenceSequence.CODING_DNA
    token_type = 'CodingDNASubstitution'
    so_id = 'SO:0001483'
    molecule_context = 'transcript'


class CodingDNASilentMutationToken(SingleNucleotideVariant):
    """SNV no change at the coding DNA reference sequence."""

    reference_sequence = ReferenceSequence.CODING_DNA
    new_nucleotide = '='
    token_type = 'CodingDNASilentMutation'
    so_id = 'SO:0002073'
    molecule_context = 'transcript'


class GenomicSubstitutionToken(SingleNucleotideVariant):
    """SNV substitution at the linear genomic reference sequence."""

    reference_sequence = ReferenceSequence.LINEAR_GENOMIC
    token_type = 'GenomicSubstitution'
    so_id = 'SO:0001483'
    molecule_context = 'genomic'


class GenomicSilentMutationToken(SingleNucleotideVariant):
    """SNV no change at the linear genomic reference sequence."""

    reference_sequence = ReferenceSequence.LINEAR_GENOMIC
    new_nucleotide = '='
    token_type = 'GenomicSilentMutation'
    so_id = 'SO:0002073'
    molecule_context = 'genomic'


class DelIns(Token):
    """A sequence alteration which included an insertion and a deletion,
    affecting 2 or more bases.
    """

    start_pos_del: Optional[str]
    end_pos_del: str
    inserted_sequence1: str
    inserted_sequence2: Optional[str]
    token_type: str
    reference_sequence: ReferenceSequence
    so_id = 'SO:1000032'
    molecule_context: str


class AminoAcidDelInsToken(Token):
    """DelIns at the protein reference sequence."""

    start_aa_del: str
    start_pos_del: int
    end_aa_del: Optional[str]
    end_pos_del: Optional[int]
    inserted_sequence: str
    reference_sequence = ReferenceSequence.PROTEIN
    so_id = 'SO:1000032'
    molecule_context = 'protein'
    token_type = 'AminoAcidDelIns'


class CodingDNADelInsToken(DelIns):
    """DelIns at the coding DNA reference sequence."""

    reference_sequence = ReferenceSequence.CODING_DNA
    token_type = 'CodingDNADelIns'
    molecule_context = 'transcript'


class GenomicDelInsToken(DelIns):
    """DelIns at the linear genomic reference sequence."""

    reference_sequence = ReferenceSequence.LINEAR_GENOMIC
    token_type = 'GenomicDelIns'
    molecule_context = 'genomic'


class LocusReferenceGenomicToken(Token):
    """Contain stable reference sequences that are used for reporting
    sequence variants with clinical implications.
    """

    id: int
    t: Optional[int]
    p: Optional[int]
    token_type = 'LocusReferenceGenomic'


class Deletion(Token):
    """The point at which one or more contiguous nucleotides were excised.
    - Sequence Ontology
    """

    start_pos_del: int
    end_pos_del: Optional[int]
    reference_sequence: ReferenceSequence
    token_type: str
    so_id: str
    molecule_context: str


class AminoAcidDeletionToken(Deletion):
    """A sequence change between the translation initiation (start) and
    termination (stop) codon where, compared to a reference sequence, one or
    more amino acids are not present (deleted) - HGVS Nomenclature
    """

    start_aa_del: str
    end_aa_del: Optional[str]
    reference_sequence = ReferenceSequence.PROTEIN
    token_type = 'AminoAcidDeletion'
    so_id = 'SO:0001604'
    molecule_context = 'protein'


class CodingDNADeletionToken(Deletion):
    """A sequence change where, compared to a reference sequence, one or
    more nucleotides are not present (deleted). - HGVS Nomenclature
    """

    reference_sequence = ReferenceSequence.CODING_DNA
    deleted_sequence: Optional[str]
    token_type = 'CodingDNADeletion'
    so_id = 'SO:0000159'
    molecule_context = 'transcript'


class GenomicDeletionToken(Deletion):
    """A sequence change where, compared to a reference sequence, one or
    more nucleotides are not present (deleted). - HGVS Nomenclature
    """

    reference_sequence = ReferenceSequence.LINEAR_GENOMIC
    deleted_sequence: Optional[str]
    token_type = 'GenomicDeletion'
    so_id = 'SO:0000159'
    molecule_context = 'genomic'
