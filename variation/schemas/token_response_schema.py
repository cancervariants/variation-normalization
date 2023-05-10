"""Module for Token Schema."""
from typing import List, Union, Dict, Any, Type, Optional, Literal
from enum import Enum

from pydantic import BaseModel
from ga4gh.vrsatile.pydantic.vrsatile_models import GeneDescriptor


class TokenType(str, Enum):
    """Define token types."""

    AMPLIFICATION = "Amplification"
    CHROMOSOME = "Chromosome"
    CODING_DNA_DELETION = "CodingDNADeletion"
    CODING_DNA_DELINS = "CodingDNADelIns"
    CODING_DNA_INSERTION = "CodingDNAInsertion"
    CODING_DNA_SILENT_MUTATION = "CodingDNASilentMutation"
    CODING_DNA_SUBSTITUTION = "CodingDNASubstitution"
    GENE = "Gene"
    GENOMIC_DELETION = "GenomicDeletion"
    GENOMIC_DELETION_RANGE = "GenomicDeletionRange"
    GENOMIC_DELINS = "GenomicDelIns"
    GENOMIC_DUPLICATION = "GenomicDuplication"
    GENOMIC_DUPLICATION_RANGE = "GenomicDuplicationRange"
    GENOMIC_INSERTION = "GenomicInsertion"
    GENOMIC_SILENT_MUTATION = "GenomicSilentMutation"
    GENOMIC_SUBSTITUTION = "GenomicSubstitution"
    GENOMIC_UNCERTAIN_DELETION = "GenomicUncertainDeletion"
    HGVS = "HGVS"
    LOCUS_REFERENCE_GENOMIC = "LocusReferenceGenomic"
    POLYPEPTIDE_TRUNCATION = "PolypeptideTruncation"
    PROTEIN_DELETION = "ProteinDeletion"
    PROTEIN_DELINS = "ProteinDelIns"
    PROTEIN_INSERTION = "ProteinInsertion"
    PROTEIN_SUBSTITUTION = "ProteinSubstitution"
    REFERENCE_SEQUENCE = "ReferenceSequence"
    SILENT_MUTATION = "SilentMutation"
    UNKNOWN = "Unknown"


class AltType(str, Enum):
    """Define alteration types."""

    AMPLIFICATION = "amplification"
    DELETION = "deletion"
    DELETION_RANGE = "deletion_range"
    UNCERTAIN_DELETION = "uncertain_deletion"
    DELINS = "delins"
    DUPLICATION = "duplication"
    DUPLICATION_RANGE = "duplication_range"
    INSERTION = "insertion"
    NONSENSE = "nonsense"
    SILENT_MUTATION = "silent_mutation"
    SUBSTITUTION = "substitution"
    UNCERTAIN_DUPLICATION = "uncertain_duplication"


AMBIGUOUS_REGIONS = {
    AltType.UNCERTAIN_DELETION,
    AltType.UNCERTAIN_DUPLICATION,
    AltType.DELETION_RANGE,
    AltType.DUPLICATION_RANGE
}


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
    FEATURE_AMPLIFICATION = "SO:0001880"


class Token(BaseModel):
    """A string from a given query."""

    token: str
    token_type: TokenType
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
                "input_string": "foo",
                "object_type": "Token",
                "nomenclature": None
            }


class GeneToken(Token):
    """Define model for gene symbol token."""

    matched_value: str
    token_type = TokenType.GENE
    gene_descriptor: Optional[GeneDescriptor]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["GeneToken"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "token": "BRAF",
                "token_type": "GeneSymbol",
                "input_string": "BRAF",
                "object_type": "Token",
                "matched_value": "BRAF"
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
    token_type: TokenType
    coordinate_type = CoordinateType.PROTEIN
    so_id: SequenceOntology
    molecule_context = "protein"
    alt_type: AltType


class PolypeptideTruncationToken(PolypeptideSequenceVariation):
    """A sequence variation of the CD that causes a truncation of the
    resulting polypeptide. (nonsense)
    """

    alt_protein = "*"
    token_type = TokenType.POLYPEPTIDE_TRUNCATION
    so_id = SequenceOntology.POLYPEPTIDE_TRUNCATION
    alt_type = AltType.NONSENSE

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

    token_type = TokenType.PROTEIN_SUBSTITUTION
    so_id = SequenceOntology.PROTEIN_SUBSTITUTION
    alt_type = AltType.SUBSTITUTION

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
                "input_string": "V600E",
                "object_type": "Token",
                "ref_protein": "V",
                "alt_protein": "E",
                "position": 600
            }


class SilentMutationToken(PolypeptideSequenceVariation):
    """A sequence variation that does not affect protein functions."""

    alt_protein = "="
    token_type = TokenType.SILENT_MUTATION
    so_id = SequenceOntology.SILENT_MUTATION
    alt_type = AltType.SILENT_MUTATION

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
                "input_string": "p.Cys188=",
                "object_type": "Token",
                "ref_protein": "Cys",
                "alt_protein": "=",
                "position": 188
            }


class TokenResponseSchema(BaseModel):
    """Define model for token response."""

    search_term: str
    tokens: List[
        Union[
            GeneToken,
            ProteinSubstitutionToken,
            PolypeptideTruncationToken,
            SilentMutationToken,
            Token
        ]
    ]

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
                        "input_string": "BRAF",
                        "object_type": "Token",
                        "matched_value": "BRAF"
                    },
                    {
                        "token": "V600E",
                        "token_type": "ProteinSubstitution",
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
    token_type: TokenType
    coordinate_type: CoordinateType
    so_id: SequenceOntology
    molecule_context: str
    alt_type: AltType


class CodingDNASubstitutionToken(SingleNucleotideVariation):
    """SNV substitution at the coding DNA reference sequence."""

    coordinate_type = CoordinateType.CODING_DNA
    token_type = TokenType.CODING_DNA_SUBSTITUTION
    so_id = SequenceOntology.SNV
    molecule_context = "transcript"
    alt_type = AltType.SUBSTITUTION


class CodingDNASilentMutationToken(SingleNucleotideVariation):
    """SNV no change at the coding DNA reference sequence."""

    coordinate_type = CoordinateType.CODING_DNA
    new_nucleotide = "="
    token_type = TokenType.CODING_DNA_SILENT_MUTATION
    so_id = SequenceOntology.NO_SEQUENCE_ALTERATION
    molecule_context = "transcript"
    alt_type = AltType.SILENT_MUTATION


class GenomicSubstitutionToken(SingleNucleotideVariation):
    """SNV substitution at the linear genomic reference sequence."""

    coordinate_type = CoordinateType.LINEAR_GENOMIC
    token_type = TokenType.GENOMIC_SUBSTITUTION
    so_id = SequenceOntology.SNV
    molecule_context = "genomic"
    alt_type = AltType.SUBSTITUTION


class GenomicSilentMutationToken(SingleNucleotideVariation):
    """SNV no change at the linear genomic reference sequence."""

    coordinate_type = CoordinateType.LINEAR_GENOMIC
    new_nucleotide = "="
    token_type = TokenType.GENOMIC_SILENT_MUTATION
    so_id = SequenceOntology.NO_SEQUENCE_ALTERATION
    molecule_context = "genomic"
    alt_type = AltType.SILENT_MUTATION


class DelIns(Token):
    """A sequence alteration which included an insertion and a deletion,
    affecting 2 or more bases.
    """

    start_pos_del: str
    end_pos_del: Optional[str]
    inserted_sequence1: str
    inserted_sequence2: Optional[str]
    token_type: TokenType
    coordinate_type: CoordinateType
    so_id = SequenceOntology.DELINS
    molecule_context: str
    alt_type = AltType.DELINS


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
    token_type = TokenType.PROTEIN_DELINS
    alt_type = AltType.DELINS


class CodingDNADelInsToken(DelIns):
    """DelIns at the coding DNA reference sequence."""

    coordinate_type = CoordinateType.CODING_DNA
    token_type = TokenType.CODING_DNA_DELINS
    molecule_context = "transcript"


class GenomicDelInsToken(DelIns):
    """DelIns at the linear genomic reference sequence."""

    coordinate_type = CoordinateType.LINEAR_GENOMIC
    token_type = TokenType.GENOMIC_DELINS
    molecule_context = "genomic"


class LocusReferenceGenomicToken(Token):
    """Contain stable reference sequences that are used for reporting
    sequence variations with clinical implications.
    """

    id: int
    t: Optional[int]
    p: Optional[int]
    token_type = TokenType.LOCUS_REFERENCE_GENOMIC


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
    token_type: TokenType
    so_id: SequenceOntology
    molecule_context: str
    alt_type = AltType.INSERTION


class ProteinInsertionToken(Insertion):
    """Protein Insertion."""

    start_aa_flank: str
    end_aa_flank: str
    coordinate_type = CoordinateType.PROTEIN
    token_type = TokenType.PROTEIN_INSERTION
    so_id = SequenceOntology.PROTEIN_INSERTION
    molecule_context = "protein"


class CodingDNAInsertionToken(Insertion):
    """Coding DNA Insertion."""

    coordinate_type = CoordinateType.CODING_DNA
    inserted_sequence2: Optional[str]
    token_type = TokenType.CODING_DNA_INSERTION
    so_id = SequenceOntology.INSERTION
    molecule_context = "transcript"


class GenomicInsertionToken(Insertion):
    """Genomic Insertion."""

    coordinate_type = CoordinateType.LINEAR_GENOMIC
    inserted_sequence2: Optional[str]
    token_type = TokenType.GENOMIC_INSERTION
    so_id = SequenceOntology.INSERTION
    molecule_context = "genomic"


class Deletion(Token):
    """The point at which one or more contiguous nucleotides were excised.
    - Sequence Ontology
    """

    start_pos_del: int
    end_pos_del: Optional[int]
    coordinate_type: CoordinateType
    token_type: TokenType
    so_id: SequenceOntology
    molecule_context: str
    alt_type: AltType.DELETION = AltType.DELETION


class ProteinDeletionToken(Deletion):
    """A sequence change between the translation initiation (start) and
    termination (stop) codon where, compared to a reference sequence, one or
    more amino acids are not present (deleted) - varnomen.hgvs.org
    """

    start_aa_del: str
    end_aa_del: Optional[str]
    deleted_aa: Optional[str]
    coordinate_type = CoordinateType.PROTEIN
    token_type = TokenType.PROTEIN_DELETION
    so_id = SequenceOntology.PROTEIN_DELETION
    molecule_context = "protein"


class CodingDNADeletionToken(Deletion):
    """A sequence change where, compared to a reference sequence, one or
    more nucleotides are not present (deleted). - varnomen.hgvs.org
    """

    coordinate_type = CoordinateType.CODING_DNA
    deleted_sequence: Optional[str]
    token_type = TokenType.CODING_DNA_DELETION
    so_id = SequenceOntology.DELETION
    molecule_context = "transcript"


class GenomicDeletionToken(Deletion):
    """A sequence change where, compared to a reference sequence, one or
    more nucleotides are not present (deleted). - varnomen.hgvs.org
    """

    coordinate_type = CoordinateType.LINEAR_GENOMIC
    deleted_sequence: Optional[str]
    token_type = TokenType.GENOMIC_DELETION
    so_id = SequenceOntology.DELETION
    molecule_context = "genomic"


class DeletionRange(Token):
    """Deletions of the form (pos_pos)_(pos_pos)."""

    start_pos1_del: Union[int, str]
    start_pos2_del: int
    end_pos1_del: int
    end_pos2_del: Union[int, str]
    token_type: TokenType
    so_id = SequenceOntology.COPY_NUMBER_LOSS
    molecule_context: str
    alt_type: Union[Literal[AltType.DELETION_RANGE], Literal[AltType.UNCERTAIN_DELETION]] = AltType.DELETION_RANGE  # noqa: E501


class GenomicDeletionRangeToken(DeletionRange):
    """Genomic deletion range token."""

    token_type = TokenType.GENOMIC_DELETION_RANGE
    molecule_context = "genomic"
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class UncertainDeletion(DeletionRange):
    """Uncertain Deletion."""

    start_pos1_del: Optional[Union[Literal["?"], int]]
    start_pos2_del: Optional[int]
    end_pos1_del: int
    end_pos2_del: Optional[Union[Literal["?"], int]]
    token_type: TokenType
    molecule_context: str
    alt_type: Literal[AltType.UNCERTAIN_DELETION] = AltType.UNCERTAIN_DELETION


class GenomicUncertainDeletionToken(UncertainDeletion):
    """Genomic uncertain deletion."""

    token_type = TokenType.GENOMIC_UNCERTAIN_DELETION
    molecule_context = "genomic"
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class Duplication(Token):
    """Duplications."""

    start_pos1_dup: Union[Literal["?"], int]
    start_pos2_dup: Optional[int]
    token_type: TokenType
    so_id = SequenceOntology.DUPLICATION
    molecule_context: str
    alt_type: Union[
        Literal[AltType.DUPLICATION],
        Literal[AltType.DUPLICATION_RANGE],
        Literal[AltType.UNCERTAIN_DUPLICATION]
    ]


class GenomicDuplicationToken(Duplication):
    """Genomic duplication token schema."""

    token_type = TokenType.GENOMIC_DUPLICATION
    molecule_context = "genomic"
    coordinate_type = CoordinateType.LINEAR_GENOMIC
    alt_type: Literal[AltType.DUPLICATION] = AltType.DUPLICATION


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
    token_type = TokenType.CHROMOSOME


class AmplificationToken(Token):
    """Amplification token"""

    token_type = TokenType.AMPLIFICATION
    molecule_context = "genomic"
    so_id = SequenceOntology.FEATURE_AMPLIFICATION
    alt_type = AltType.AMPLIFICATION
