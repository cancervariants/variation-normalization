"""Module for schemas related to tokenization."""
from enum import Enum
from typing import Optional

from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.vrsatile.pydantic.vrsatile_models import GeneDescriptor
from pydantic import BaseModel, StrictInt, StrictStr

from variation.schemas.app_schemas import AmbiguousRegexType
from variation.schemas.variation_schema import (
    Deletion,
    DelIns,
    DupDelAmbiguous,
    Duplication,
    Insertion,
    ProteinDeletion,
    ProteinDelIns,
    ProteinInsertion,
    ProteinReferenceAgree,
    ReferenceAgree,
    StopGain,
    Substitution,
)


class TokenType(str, Enum):
    """Define token types."""

    AMPLIFICATION = "amplification"
    CDNA_DELETION = "cdna_deletion"
    CDNA_DELINS = "cdna_delins"
    CDNA_INSERTION = "cdna_insertion"
    CDNA_REFERENCE_AGREE = "cdna_reference_agree"
    CDNA_SUBSTITUTION = "cdna_substitution"
    GENE = "gene"
    GENOMIC_DELETION = "genomic_deletion"
    GENOMIC_DELETION_AMBIGUOUS = "genomic_deletion_ambiguous"
    GENOMIC_DELINS = "genomic_delins"
    GENOMIC_DUPLICATION = "genomic_duplication"
    GENOMIC_DUPLICATION_AMBIGUOUS = "genomic_duplication_ambiguous"
    GENOMIC_INSERTION = "genomic_insertion"
    GENOMIC_REFERENCE_AGREE = "genomic_reference_agree"
    GENOMIC_SUBSTITUTION = "genomic_substitution"
    GNOMAD_VCF = "gnomad_vcf"
    HGVS = "hgvs"
    PROTEIN_STOP_GAIN = "protein_stop_gain"
    PROTEIN_DELETION = "protein_deletion"
    PROTEIN_DELINS = "protein_delins"
    PROTEIN_INSERTION = "protein_insertion"
    PROTEIN_SUBSTITUTION = "protein_substitution"
    PROTEIN_REFERENCE_AGREE = "protein_reference_agree"
    UNKNOWN = "unknown"


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
AMBIGUOUS_REGIONS = {AltType.DELETION_AMBIGUOUS, AltType.DUPLICATION_AMBIGUOUS}


class Token(BaseModel):
    """A string from a given query."""

    token: StrictStr
    token_type: TokenType
    input_string: StrictStr


class HgvsToken(Token):
    """HGVS Token"""

    token_type = TokenType.HGVS
    accession: StrictStr
    coordinate_type: AnnotationLayer
    change: StrictStr


class GnomadVcfToken(Token):
    """Gnomad VCF Token"""

    token_type = TokenType.GNOMAD_VCF
    coordinate_type = AnnotationLayer.GENOMIC
    chromosome: StrictStr
    pos: StrictInt
    ref: StrictStr
    alt: StrictStr


class GenomicSubstitutionToken(Token, Substitution):
    """Genomic substitution token"""

    token_type = TokenType.GENOMIC_SUBSTITUTION
    coordinate_type = AnnotationLayer.GENOMIC


class CdnaSubstitutionToken(Token, Substitution):
    """Token for substitution on cDNA reference sequence"""

    token_type = TokenType.CDNA_SUBSTITUTION
    coordinate_type = AnnotationLayer.CDNA


class ProteinSubstitutionToken(Token, Substitution):
    """Token for substitution on protein reference sequence"""

    token_type = TokenType.PROTEIN_SUBSTITUTION
    coordinate_type = AnnotationLayer.PROTEIN


class ProteinStopGainToken(Token, StopGain):
    """Token for stop gain on protein reference sequence"""

    token_type = TokenType.PROTEIN_STOP_GAIN
    coordinate_type = AnnotationLayer.PROTEIN


class ProteinReferenceAgreeToken(Token, ProteinReferenceAgree):
    """Token for reference agree on protein reference sequence"""

    coordinate_type = AnnotationLayer.PROTEIN
    token_type = TokenType.PROTEIN_REFERENCE_AGREE


class CdnaReferenceAgreeToken(Token, ReferenceAgree):
    """Token for reference agree on cDNA reference sequence"""

    coordinate_type = AnnotationLayer.CDNA
    token_type = TokenType.CDNA_REFERENCE_AGREE


class GenomicReferenceAgreeToken(Token, ReferenceAgree):
    """Token for reference agree on genomic reference sequence"""

    coordinate_type = AnnotationLayer.GENOMIC
    token_type = TokenType.GENOMIC_REFERENCE_AGREE


class ProteinDeletionToken(Token, ProteinDeletion):
    """Token for deletion on protein reference sequence"""

    token_type = TokenType.PROTEIN_DELETION
    coordinate_type = AnnotationLayer.PROTEIN


class CdnaDeletionToken(Token, Deletion):
    """Token for deletion on cdna reference sequence"""

    token_type = TokenType.CDNA_DELETION
    coordinate_type = AnnotationLayer.CDNA


class GenomicDeletionToken(Token, Deletion):
    """Token for deletion on genomic reference sequence"""

    token_type = TokenType.GENOMIC_DELETION
    coordinate_type = AnnotationLayer.GENOMIC


class GenomicDeletionAmbiguousToken(Token, DupDelAmbiguous):
    """Token for ambiguous deletion on genomic reference sequence"""

    token_type = TokenType.GENOMIC_DELETION_AMBIGUOUS
    coordinate_type = AnnotationLayer.GENOMIC
    ambiguous_regex_type: AmbiguousRegexType


class ProteinDelInsToken(Token, ProteinDelIns):
    """Token for delins on protein reference sequence"""

    token_type = TokenType.PROTEIN_DELINS
    coordinate_type = AnnotationLayer.PROTEIN


class CdnaDelInsToken(Token, DelIns):
    """Token for delins on cdna reference sequence"""

    token_type = TokenType.CDNA_DELINS
    coordinate_type = AnnotationLayer.CDNA


class GenomicDelInsToken(Token, DelIns):
    """Token for delins on genomic reference sequence"""

    token_type = TokenType.GENOMIC_DELINS
    coordinate_type = AnnotationLayer.GENOMIC


class CdnaInsertionToken(Token, Insertion):
    """Token for insertion on cdna reference sequence"""

    token_type = TokenType.CDNA_INSERTION
    coordinate_type = AnnotationLayer.CDNA


class GenomicInsertionToken(Token, Insertion):
    """Token for insertion on genomic reference sequence"""

    token_type = TokenType.GENOMIC_INSERTION
    coordinate_type = AnnotationLayer.GENOMIC


class ProteinInsertionToken(Token, ProteinInsertion):
    """Token for insertion on protein reference sequence"""

    token_type = TokenType.PROTEIN_INSERTION
    coordinate_type = AnnotationLayer.PROTEIN


class GenomicDuplicationToken(Token, Duplication):
    """Duplication on genomic reference sequence"""

    token_type = TokenType.GENOMIC_DUPLICATION
    coordinate_type = AnnotationLayer.GENOMIC


class GenomicDuplicationAmbiguousToken(Token, DupDelAmbiguous):
    """Ambiguous duplication on genomic reference sequence"""

    token_type = TokenType.GENOMIC_DUPLICATION_AMBIGUOUS
    coordinate_type = AnnotationLayer.GENOMIC
    ambiguous_regex_type: AmbiguousRegexType


class AmplificationToken(Token):
    """Token for amplification"""

    token_type = TokenType.AMPLIFICATION


class GeneToken(Token):
    """Token for genes"""

    matched_value: StrictStr
    token_type = TokenType.GENE
    gene_descriptor: Optional[GeneDescriptor]
