"""Module for Token Schema."""
from typing import List, Union, Dict, Any, Type, Optional, Literal
from enum import IntEnum, Enum

from pydantic import BaseModel


class TokenMatchType(IntEnum):
    """Gene Token match types."""

    ID = 1
    SYMBOL = 2
    PREVIOUS = 3
    ALIAS = 4
    UNSPECIFIED = 5


class TokenType(str, Enum):
    """Define token types."""

    # TODO: Add other token types (issue 162)
    GENOMIC_DUPLICATION = "GenomicDuplication"
    GENOMIC_DUPLICATION_RANGE = "GenomicDuplicationRange"


class Nomenclature(str, Enum):
    """Define nomenclatures that are supported"""

    FREE_TEXT = "free_text"
    HGVS = "hgvs"
    GNOMAD_VCF = "gnomad_vcf"


class SequenceOntology(str, Enum):
    """Define SequenceOntology codes"""

    POLYPEPTIDE_TRUNCATION = "SO:0001617"
    PROTEIN_SUBSTITUTION = "SO:0001606"
    SILENT_MUTATION = "SO:0001017"
    SNV = "SO:0001483"
    NO_SEQUENCE_ALTERATION = "SO:0002073"
    DELINS = "SO:1000032"
    PROTEIN_INSERTION = "SO:0001605"
    INSERTION = "SO:0000667"
    PROTEIN_DELETION = "SO:0001604"
    DELETION = "SO:0000159"
    COPY_NUMBER_LOSS = "SO:0001743"
    DUPLICATION = "SO:1000035"
    COPY_NUMBER_GAIN = "SO:0001742"


class Token(BaseModel):
    """A string from a given query."""

    token: str
    token_type: str
    match_type: TokenMatchType
    input_string: str
    object_type = "Token"
    nomenclature: Optional[Nomenclature]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["Token"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "token": "",
                "token_type": "Unknown",
                "match_type": 5,
                "input_string": "foo",
                "object_type": "Token",
                "nomenclature": None
            }


class GeneMatchToken(Token):
    """Define model for gene symbol token."""

    matched_value: str
    token_type = "GeneSymbol"

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["GeneMatchToken"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
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
    token_type = "GenePair"

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["GenePairMatchToken"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
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


class CoordinateType(str, Enum):
    """Define constraints for coordinate types."""

    CODING_DNA = "c"
    LINEAR_GENOMIC = "g"
    PROTEIN = "p"


class PolypeptideSequenceVariation(Token):
    """Polypeptide Sequence Variation Token Class."""

    ref_protein: str
    alt_protein: str
    position: int
    token_type: str
    coordinate_type = CoordinateType.PROTEIN
    so_id: SequenceOntology
    molecule_context = "protein"
    alt_type: str


class PolypeptideTruncationToken(PolypeptideSequenceVariation):
    """A sequence variation of the CD that causes a truncation of the
    resulting polypeptide. (nonsense)
    """

    alt_protein = "*"
    token_type = "PolypeptideTruncation"
    so_id = SequenceOntology.POLYPEPTIDE_TRUNCATION
    alt_type = "nonsense"

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["PolypeptideTruncationToken"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "token": "Tyr365Ter",
                "token_type": "PolypeptideTruncation",
                "match_type": 5,
                "input_string": "Tyr365Ter",
                "object_type": "Token",
                "ref_protein": "Tyr",
                "alt_protein": "Ter",
                "position": 365
            }


class ProteinSubstitutionToken(PolypeptideSequenceVariation):
    """A sequence variation of a codon resulting in the substitution of one
    amino acid for another in the resulting polypeptide. (missense)
    """

    token_type = "ProteinSubstitution"
    so_id = SequenceOntology.PROTEIN_SUBSTITUTION
    alt_type = "substitution"

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["ProteinSubstitutionToken"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "token": "V600E",
                "token_type": "ProteinSubstitution",
                "match_type": 5,
                "input_string": "V600E",
                "object_type": "Token",
                "ref_protein": "V",
                "alt_protein": "E",
                "position": 600
            }


class SilentMutationToken(PolypeptideSequenceVariation):
    """A sequence variation that does not affect protein functions."""

    alt_protein = "="
    token_type = "SilentMutation"
    so_id = SequenceOntology.SILENT_MUTATION
    alt_type = "silent_mutation"

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["SilentMutationToken"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
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
                       ProteinSubstitutionToken, PolypeptideTruncationToken,
                       SilentMutationToken, Token]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["TokenResponseSchema"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
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
                        "token_type": "ProteinSubstitution",
                        "match_type": 5,
                        "input_string": "V600E",
                        "object_type": "Token",
                        "ref_protein": "V",
                        "alt_protein": "E",
                        "position": 600
                    }
                ]
            }


class SingleNucleotideVariation(Token):
    """Single nucleotide positions in genomic DNA at which different
    sequence alternatives exist.
    """

    position: int
    ref_nucleotide: Optional[str]
    new_nucleotide: str
    token_type: str
    coordinate_type: CoordinateType
    so_id: SequenceOntology
    molecule_context: str
    alt_type: str


class CodingDNASubstitutionToken(SingleNucleotideVariation):
    """SNV substitution at the coding DNA reference sequence."""

    coordinate_type = CoordinateType.CODING_DNA
    token_type = "CodingDNASubstitution"
    so_id = SequenceOntology.SNV
    molecule_context = "transcript"
    alt_type = "substitution"


class CodingDNASilentMutationToken(SingleNucleotideVariation):
    """SNV no change at the coding DNA reference sequence."""

    coordinate_type = CoordinateType.CODING_DNA
    new_nucleotide = "="
    token_type = "CodingDNASilentMutation"
    so_id = SequenceOntology.NO_SEQUENCE_ALTERATION
    molecule_context = "transcript"
    alt_type = "silent_mutation"


class GenomicSubstitutionToken(SingleNucleotideVariation):
    """SNV substitution at the linear genomic reference sequence."""

    coordinate_type = CoordinateType.LINEAR_GENOMIC
    token_type = "GenomicSubstitution"
    so_id = SequenceOntology.SNV
    molecule_context = "genomic"
    alt_type = "substitution"


class GenomicSilentMutationToken(SingleNucleotideVariation):
    """SNV no change at the linear genomic reference sequence."""

    coordinate_type = CoordinateType.LINEAR_GENOMIC
    new_nucleotide = "="
    token_type = "GenomicSilentMutation"
    so_id = SequenceOntology.NO_SEQUENCE_ALTERATION
    molecule_context = "genomic"
    alt_type = "silent_mutation"


class DelIns(Token):
    """A sequence alteration which included an insertion and a deletion,
    affecting 2 or more bases.
    """

    start_pos_del: str
    end_pos_del: Optional[str]
    inserted_sequence1: str
    inserted_sequence2: Optional[str]
    token_type: str
    coordinate_type: CoordinateType
    so_id = SequenceOntology.DELINS
    molecule_context: str
    alt_type = "delins"


class ProteinDelInsToken(Token):
    """DelIns at the protein reference sequence."""

    start_aa_del: str
    start_pos_del: int
    end_aa_del: Optional[str]
    end_pos_del: Optional[int]
    inserted_sequence: str
    coordinate_type = CoordinateType.PROTEIN
    so_id = SequenceOntology.DELINS
    molecule_context = "protein"
    token_type = "ProteinDelIns"
    alt_type = "delins"


class CodingDNADelInsToken(DelIns):
    """DelIns at the coding DNA reference sequence."""

    coordinate_type = CoordinateType.CODING_DNA
    token_type = "CodingDNADelIns"
    molecule_context = "transcript"


class GenomicDelInsToken(DelIns):
    """DelIns at the linear genomic reference sequence."""

    coordinate_type = CoordinateType.LINEAR_GENOMIC
    token_type = "GenomicDelIns"
    molecule_context = "genomic"


class LocusReferenceGenomicToken(Token):
    """Contain stable reference sequences that are used for reporting
    sequence variations with clinical implications.
    """

    id: int
    t: Optional[int]
    p: Optional[int]
    token_type = "LocusReferenceGenomic"


class Insertion(Token):
    """a sequence change between the translation initiation (start) and
    termination (stop) codon where, compared to the reference sequence,
    one or more amino acids are inserted, which is not a frame shift and
    where the insertion is not a copy of a sequence immediately N-terminal
    (5') - varnomen.hgvs.org
    """

    start_pos_flank: int
    end_pos_flank: int
    inserted_sequence: str
    coordinate_type: CoordinateType
    token_type: str
    so_id: SequenceOntology
    molecule_context: str
    alt_type = "insertion"


class ProteinInsertionToken(Insertion):
    """Protein Insertion."""

    start_aa_flank: str
    end_aa_flank: str
    coordinate_type = CoordinateType.PROTEIN
    token_type = "ProteinInsertion"
    so_id = SequenceOntology.PROTEIN_INSERTION
    molecule_context = "protein"


class CodingDNAInsertionToken(Insertion):
    """Coding DNA Insertion."""

    coordinate_type = CoordinateType.CODING_DNA
    inserted_sequence2: Optional[str]
    token_type = "CodingDNAInsertion"
    so_id = SequenceOntology.INSERTION
    molecule_context = "transcript"


class GenomicInsertionToken(Insertion):
    """Genomic Insertion."""

    coordinate_type = CoordinateType.LINEAR_GENOMIC
    inserted_sequence2: Optional[str]
    token_type = "GenomicInsertion"
    so_id = SequenceOntology.INSERTION
    molecule_context = "genomic"


class DeletionAltType(str, Enum):
    """Define alt types for deletions."""

    DELETION = "deletion"
    DELETION_RANGE = "deletion_range"
    UNCERTAIN_DELETION = "uncertain_deletion"


class Deletion(Token):
    """The point at which one or more contiguous nucleotides were excised.
    - Sequence Ontology
    """

    start_pos_del: int
    end_pos_del: Optional[int]
    coordinate_type: CoordinateType
    token_type: str
    so_id: SequenceOntology
    molecule_context: str
    alt_type: DeletionAltType.DELETION = DeletionAltType.DELETION


class ProteinDeletionToken(Deletion):
    """A sequence change between the translation initiation (start) and
    termination (stop) codon where, compared to a reference sequence, one or
    more amino acids are not present (deleted) - varnomen.hgvs.org
    """

    start_aa_del: str
    end_aa_del: Optional[str]
    coordinate_type = CoordinateType.PROTEIN
    token_type = "ProteinDeletion"
    so_id = SequenceOntology.PROTEIN_DELETION
    molecule_context = "protein"


class CodingDNADeletionToken(Deletion):
    """A sequence change where, compared to a reference sequence, one or
    more nucleotides are not present (deleted). - varnomen.hgvs.org
    """

    coordinate_type = CoordinateType.CODING_DNA
    deleted_sequence: Optional[str]
    token_type = "CodingDNADeletion"
    so_id = SequenceOntology.DELETION
    molecule_context = "transcript"


class GenomicDeletionToken(Deletion):
    """A sequence change where, compared to a reference sequence, one or
    more nucleotides are not present (deleted). - varnomen.hgvs.org
    """

    coordinate_type = CoordinateType.LINEAR_GENOMIC
    deleted_sequence: Optional[str]
    token_type = "GenomicDeletion"
    so_id = SequenceOntology.DELETION
    molecule_context = "genomic"


class DeletionRange(Token):
    """Deletions of the form (pos_pos)_(pos_pos)."""

    start_pos1_del: Union[int, str]
    start_pos2_del: int
    end_pos1_del: int
    end_pos2_del: Union[int, str]
    token_type: str
    so_id = SequenceOntology.COPY_NUMBER_LOSS
    molecule_context: str
    alt_type: Union[Literal[DeletionAltType.DELETION_RANGE], Literal[DeletionAltType.UNCERTAIN_DELETION]] = DeletionAltType.DELETION_RANGE  # noqa: E501


class GenomicDeletionRangeToken(DeletionRange):
    """Genomic deletion range token."""

    token_type = "GenomicDeletionRange"
    molecule_context = "genomic"
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class UncertainDeletion(DeletionRange):
    """Uncertain Deletion."""

    start_pos1_del: Optional[Union[Literal["?"], int]]
    start_pos2_del: Optional[int]
    end_pos1_del: int
    end_pos2_del: Optional[Union[Literal["?"], int]]
    token_type: str
    molecule_context: str
    alt_type: Literal[DeletionAltType.UNCERTAIN_DELETION] = DeletionAltType.UNCERTAIN_DELETION  # noqa: E501


class GenomicUncertainDeletionToken(UncertainDeletion):
    """Genomic uncertain deletion."""

    token_type = "GenomicUncertainDeletion"
    molecule_context = "genomic"
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class DuplicationAltType(str, Enum):
    """Define alt types for duplications."""

    DUPLICATION = "duplication"
    DUPLICATION_RANGE = "duplication_range"
    UNCERTAIN_DUPLICATION = "uncertain_duplication"


class Duplication(Token):
    """Duplications."""

    start_pos1_dup: Union[Literal["?"], int]
    start_pos2_dup: Optional[int]
    token_type: TokenType
    so_id = SequenceOntology.DUPLICATION
    molecule_context: str
    alt_type: DuplicationAltType


class GenomicDuplicationToken(Duplication):
    """Genomic duplication token schema."""

    token_type = TokenType.GENOMIC_DUPLICATION
    molecule_context = "genomic"
    coordinate_type = CoordinateType.LINEAR_GENOMIC
    alt_type: Literal[DuplicationAltType.DUPLICATION] = DuplicationAltType.DUPLICATION  # noqa: E501


class DuplicationRange(Duplication):
    """Duplications of the form (#_#)_(#_#)dup"""

    end_pos1_dup: int
    end_pos2_dup: Optional[Union[Literal["?"], int]]
    so_id = SequenceOntology.COPY_NUMBER_GAIN


class GenomicDuplicationRangeToken(DuplicationRange):
    """Genomic Duplication Range token schema"""

    token_type = TokenType.GENOMIC_DUPLICATION_RANGE
    molecule_context = "genomic"
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class ChromosomeToken(Token):
    """Chromosome token"""

    chromosome: str
    token_type = "Chromosome"
