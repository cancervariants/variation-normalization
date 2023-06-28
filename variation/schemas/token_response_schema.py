"""Module for schemas related to tokenization."""
from typing import List, Optional
from enum import Enum

from pydantic import BaseModel, StrictInt, StrictStr
from ga4gh.vrsatile.pydantic.vrsatile_models import GeneDescriptor

from variation.schemas.variation_schema import (
    DelIns, Deletion, ProteinDelIns, ProteinInsertion, Substitution, Insertion,
    ProteinDeletion, ReferenceAgree, ProteinReferenceAgree, StopGain, Duplication,
    DupDelAmbiguous
)


class TokenType(str, Enum):
    """Define token types."""

    AMPLIFICATION = "Amplification"
    CODING_DNA_DELETION = "CodingDNADeletion"
    CODING_DNA_DELINS = "CodingDNADelIns"
    CODING_DNA_INSERTION = "CodingDNAInsertion"
    CODING_DNA_REFERENCE_AGREE = "CodingDNAReferenceAgree"
    CODING_DNA_SUBSTITUTION = "CodingDNASubstitution"
    GENE = "Gene"
    GENOMIC_DELETION = "GenomicDeletion"
    GENOMIC_DELETION_AMBIGUOUS = "GenomicDeletionAmbiguous"
    GENOMIC_DELINS = "GenomicDelIns"
    GENOMIC_DUPLICATION = "GenomicDuplication"
    GENOMIC_DUPLICATION_AMBIGUOUS = "GenomicDuplicationAmbiguous"
    GENOMIC_INSERTION = "GenomicInsertion"
    GENOMIC_REFERENCE_AGREE = "GenomicReferenceAgree"
    GENOMIC_SUBSTITUTION = "GenomicSubstitution"
    GNOMAD_VCF = "GnomadVcf"
    HGVS = "Hgvs"
    HGVS_SUBSTITUTION = "HgvsSubstitution"
    PROTEIN_STOP_GAIN = "ProteinStopGain"
    PROTEIN_DELETION = "ProteinDeletion"
    PROTEIN_DELINS = "ProteinDelIns"
    PROTEIN_INSERTION = "ProteinInsertion"
    PROTEIN_SUBSTITUTION = "ProteinSubstitution"
    PROTEIN_REFERENCE_AGREE = "ProteinReferenceAgree"
    UNKNOWN = "Unknown"


class AltType(str, Enum):
    """Define alteration types."""

    AMPLIFICATION = "amplification"
    DELETION = "deletion"
    DELETION_AMBIGUOUS = "deletion_ambiguous"
    DELINS = "delins"
    DUPLICATION = "duplication"
    DUPLICATION_AMBIGUOUS = "duplication_ambiguous"
    INSERTION = "insertion"
    NONSENSE = "nonsense"
    REFERENCE_AGREE = "reference_agree"
    SUBSTITUTION = "substitution"
    STOP_GAIN = "stop_gain"


# Ambiguous region alt types
AMBIGUOUS_REGIONS = {
    AltType.DELETION_AMBIGUOUS,
    AltType.DUPLICATION_AMBIGUOUS
}


class CoordinateType(str, Enum):
    """Define constraints for coordinate types."""

    CODING_DNA = "c"
    LINEAR_GENOMIC = "g"
    PROTEIN = "p"


class Token(BaseModel):
    """A string from a given query."""

    token: str
    token_type: TokenType
    input_string: str


class HgvsToken(Token):
    """HGVS Token"""

    token_type = TokenType.HGVS
    accession: StrictStr
    coordinate_type: CoordinateType
    change: StrictStr


class GnomadVcfToken(Token):
    """Gnomad VCF Token"""

    token_type = TokenType.GNOMAD_VCF
    coordinate_type = CoordinateType.LINEAR_GENOMIC
    chromosome: StrictStr
    pos: StrictInt
    ref: StrictStr
    alt: StrictStr


class GenomicSubstitutionToken(Token, Substitution):
    """Genomic substitution token"""

    token_type = TokenType.GENOMIC_SUBSTITUTION
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class CodingDNASubstitutionToken(Token, Substitution):
    """Token for substitution on coding DNA reference sequence"""

    token_type = TokenType.CODING_DNA_SUBSTITUTION
    coordinate_type = CoordinateType.CODING_DNA


class ProteinSubstitutionToken(Token, Substitution):
    """Token for substitution on protein reference sequence"""

    token_type = TokenType.PROTEIN_SUBSTITUTION
    coordinate_type = CoordinateType.PROTEIN


class ProteinStopGainToken(Token, StopGain):
    """Token for stop gain on protein reference sequence"""

    token_type = TokenType.PROTEIN_STOP_GAIN
    coordinate_type = CoordinateType.PROTEIN


class ProteinReferenceAgreeToken(Token, ProteinReferenceAgree):
    """Token for reference agree on protein reference sequence"""

    coordinate_type = CoordinateType.PROTEIN
    token_type = TokenType.PROTEIN_REFERENCE_AGREE


class CodingDNAReferenceAgreeToken(Token, ReferenceAgree):
    """Token for reference agree on coding DNA reference sequence"""

    coordinate_type = CoordinateType.CODING_DNA
    token_type = TokenType.CODING_DNA_REFERENCE_AGREE


class GenomicReferenceAgreeToken(Token, ReferenceAgree):
    """Token for reference agree on genomic reference sequence"""

    coordinate_type = CoordinateType.LINEAR_GENOMIC
    token_type = TokenType.GENOMIC_REFERENCE_AGREE


class ProteinDeletionToken(Token, ProteinDeletion):
    """Token for deletion on protein reference sequence"""

    token_type = TokenType.PROTEIN_DELETION
    coordinate_type = CoordinateType.PROTEIN


class CodingDNADeletionToken(Token, Deletion):
    """Token for deletion on coding dna reference sequence"""

    token_type = TokenType.CODING_DNA_DELETION
    coordinate_type = CoordinateType.CODING_DNA


class GenomicDeletionToken(Token, Deletion):
    """Token for deletion on genomic reference sequence"""

    token_type = TokenType.GENOMIC_DELETION
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class GenomicDeletionAmbiguousToken(Token, DupDelAmbiguous):
    """Token for ambiguous deletion on genomic reference sequence"""

    token_type = TokenType.GENOMIC_DELETION_AMBIGUOUS
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class ProteinDelInsToken(Token, ProteinDelIns):
    """Token for delins on protein reference sequence"""

    token_type = TokenType.PROTEIN_DELINS
    coordinate_type = CoordinateType.PROTEIN


class CodingDNADelInsToken(Token, DelIns):
    """Token for delins on coding dna reference sequence"""

    token_type = TokenType.CODING_DNA_DELINS
    coordinate_type = CoordinateType.CODING_DNA


class GenomicDelInsToken(Token, DelIns):
    """Token for delins on genomic reference sequence"""

    token_type = TokenType.GENOMIC_DELINS
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class CodingDNAInsertionToken(Token, Insertion):
    """Token for insertion on coding dna reference sequence"""

    token_type = TokenType.CODING_DNA_INSERTION
    coordinate_type = CoordinateType.CODING_DNA


class GenomicInsertionToken(Token, Insertion):
    """Token for insertion on genomic reference sequence"""

    token_type = TokenType.GENOMIC_INSERTION
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class ProteinInsertionToken(Token, ProteinInsertion):
    """Token for insertion on protein reference sequence"""

    token_type = TokenType.PROTEIN_INSERTION
    coordinate_type = CoordinateType.PROTEIN


class GenomicDuplicationToken(Token, Duplication):
    """Duplication on genomic reference sequence"""

    token_type = TokenType.GENOMIC_DUPLICATION
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class GenomicDuplicationAmbiguousToken(Token, DupDelAmbiguous):
    """Ambiguous duplication on genomic reference sequence"""

    token_type = TokenType.GENOMIC_DUPLICATION_AMBIGUOUS
    coordinate_type = CoordinateType.LINEAR_GENOMIC


class AmplificationToken(Token):
    """Token for amplification"""

    token_type = TokenType.AMPLIFICATION
    alt_type = AltType.AMPLIFICATION


class GeneToken(Token):
    """Token for genes"""

    matched_value: str
    token_type = TokenType.GENE
    gene_descriptor: Optional[GeneDescriptor]


class TokenResponseSchema(BaseModel):
    """Define model for token response."""

    search_term: str
    tokens: List[Token]
