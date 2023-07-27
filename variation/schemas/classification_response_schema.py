"""Module for Classification schema."""
from typing import List, Optional, Union, Literal
from enum import Enum

from pydantic import BaseModel, StrictStr
from ga4gh.vrsatile.pydantic.vrsatile_models import MoleculeContext

from variation.schemas.token_response_schema import Token, GeneToken
from variation.schemas.variation_schema import (
    ProteinDelIns, Substitution, Deletion, Insertion, ProteinDeletion, ProteinInsertion,
    ReferenceAgree, ProteinReferenceAgree, DelIns, StopGain, Duplication,
    DupDelAmbiguous
)


class SequenceOntology(str, Enum):
    """Define SequenceOntology codes"""

    STOP_GAIN = "SO:0001587"
    PROTEIN_SUBSTITUTION = "SO:0001606"
    SNV = "SO:0001483"
    MNV = "SO:0002007"
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
    SUBSTITUTION = "SO:1000002"


class Nomenclature(str, Enum):
    """Define nomenclatures that are supported"""

    FREE_TEXT = "free_text"
    HGVS = "hgvs"
    GNOMAD_VCF = "gnomad_vcf"


class ClassificationType(str, Enum):
    """Enums for Classification Types."""

    PROTEIN_SUBSTITUTION = "protein_substitution"
    PROTEIN_STOP_GAIN = "protein_stop_gain"
    PROTEIN_REFERENCE_AGREE = "protein_reference_agree"
    PROTEIN_DELINS = "protein_delins"
    CDNA_SUBSTITUTION = "cdna_substitution"
    GENOMIC_SUBSTITUTION = "genomic_substitution"
    CDNA_REFERENCE_AGREE = "cdna_reference_agree"
    GENOMIC_REFERENCE_AGREE = "genomic_reference_agree"
    CDNA_DELINS = "cdna_delins"
    GENOMIC_DELINS = "genomic_delins"
    PROTEIN_DELETION = "protein_deletion"
    CDNA_DELETION = "cdna_deletion"
    GENOMIC_DELETION = "genomic_deletion"
    GENOMIC_DELETION_AMBIGUOUS = "genomic_deletion_ambiguous"
    PROTEIN_INSERTION = "protein_insertion"
    CDNA_INSERTION = "cdna_insertion"
    GENOMIC_INSERTION = "genomic_insertion"
    GENOMIC_DUPLICATION = "genomic_duplication"
    GENOMIC_DUPLICATION_AMBIGUOUS = "genomic_duplication_ambiguous"
    AMPLIFICATION = "amplification"


class Classification(BaseModel):
    """Classification for a list of tokens."""

    classification_type: ClassificationType
    matching_tokens: List[Token]
    nomenclature: Nomenclature
    molecule_context: MoleculeContext
    so_id: SequenceOntology
    gene_token: Optional[GeneToken]
    ac: Optional[StrictStr]


class ProteinSubstitutionClassification(Classification, Substitution):
    """Define protein substitution classification"""

    classification_type = ClassificationType.PROTEIN_SUBSTITUTION
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.PROTEIN_SUBSTITUTION


class GenomicSubstitutionClassification(Classification, Substitution):
    """Define genomic substitution classification"""

    classification_type = ClassificationType.GENOMIC_SUBSTITUTION
    molecule_context = MoleculeContext.GENOMIC
    so_id: Union[Literal[SequenceOntology.SNV], Literal[SequenceOntology.MNV]]


class CdnaSubstitutionClassification(Classification, Substitution):
    """Define cdna substitution classification"""

    classification_type = ClassificationType.CDNA_SUBSTITUTION
    molecule_context = MoleculeContext.TRANSCRIPT
    so_id: Union[Literal[SequenceOntology.SNV], Literal[SequenceOntology.MNV]]


class ProteinStopGainClassification(Classification, StopGain):
    """Define protein stop gain classification"""

    classification_type = ClassificationType.PROTEIN_STOP_GAIN
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.STOP_GAIN


class ProteinReferenceAgreeClassification(Classification, ProteinReferenceAgree):
    """Define protein reference agree classification"""

    classification_type = ClassificationType.PROTEIN_REFERENCE_AGREE
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.NO_SEQUENCE_ALTERATION


class CdnaReferenceAgreeClassification(Classification, ReferenceAgree):
    """Define cdna reference agree classification"""

    classification_type = ClassificationType.CDNA_REFERENCE_AGREE
    molecule_context = MoleculeContext.TRANSCRIPT
    so_id = SequenceOntology.NO_SEQUENCE_ALTERATION


class GenomicReferenceAgreeClassification(Classification, ReferenceAgree):
    """Define genomic reference agree classification"""

    classification_type = ClassificationType.GENOMIC_REFERENCE_AGREE
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.NO_SEQUENCE_ALTERATION


class ProteinInsertionClassification(Classification, ProteinInsertion):
    """Define protein insertion classification"""

    classification_type = ClassificationType.PROTEIN_INSERTION
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.PROTEIN_INSERTION


class CdnaInsertionClassification(Classification, Insertion):
    """Define cdna insertion classification"""

    classification_type = ClassificationType.CDNA_INSERTION
    molecule_context = MoleculeContext.TRANSCRIPT
    so_id = SequenceOntology.INSERTION


class GenomicInsertionClassification(Classification, Insertion):
    """Define genomic insertion classification"""

    classification_type = ClassificationType.GENOMIC_INSERTION
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.INSERTION


class ProteinDeletionClassification(Classification, ProteinDeletion):
    """Define protein deletion classification"""

    classification_type = ClassificationType.PROTEIN_DELETION
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.PROTEIN_DELETION


class GenomicDeletionClassification(Classification, Deletion):
    """Define genomic deletion classification"""

    classification_type = ClassificationType.GENOMIC_DELETION
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.DELETION


class CdnaDeletionClassification(Classification, Deletion):
    """Define cdna classification"""

    classification_type = ClassificationType.CDNA_DELETION
    molecule_context = MoleculeContext.TRANSCRIPT
    so_id = SequenceOntology.DELETION


class ProteinDelInsClassification(Classification, ProteinDelIns):
    """Define protein delins classification"""

    classification_type = ClassificationType.PROTEIN_DELINS
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.DELINS


class CdnaDelInsClassification(Classification, DelIns):
    """Define cdna delins classification"""

    classification_type = ClassificationType.CDNA_DELINS
    molecule_context = MoleculeContext.TRANSCRIPT
    so_id = SequenceOntology.DELINS


class GenomicDelInsClassification(Classification, DelIns):
    """Define genomic delins classification"""

    classification_type = ClassificationType.GENOMIC_DELINS
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.DELINS


class GenomicDuplicationClassification(Classification, Duplication):
    """Define genomic duplication classification"""

    classification_type = ClassificationType.GENOMIC_DUPLICATION
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.DUPLICATION


class AmbiguousType(str, Enum):
    """Define ambiguous type  which helps determine the ambiguous expression format"""

    AMBIGUOUS_1 = "(#_#)_(#_#)"
    AMBIGUOUS_2 = "(?_#)_(#_?)"
    AMBIGUOUS_3 = "(#_?)_(?_#)"  # Not yet supported
    AMBIGUOUS_4 = "(#_#)_#"  # Not yet supported
    AMBIGUOUS_5 = "(?_#)_#"
    AMBIGUOUS_6 = "#_(#_#)"  # Not yet supported
    AMBIGUOUS_7 = "#_(#_?)"


class GenomicDuplicationAmbiguousClassification(Classification, DupDelAmbiguous):
    """Define genomic duplication ambiguous classification"""

    classification_type = ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.DUPLICATION
    ambiguous_type: AmbiguousType


class GenomicDeletionAmbiguousClassification(Classification, DupDelAmbiguous):
    """Define genomic deletion ambiguous classification"""

    classification_type = ClassificationType.GENOMIC_DELETION_AMBIGUOUS
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.DELETION
    ambiguous_type: AmbiguousType


class AmplificationClassification(Classification):
    """Define amplification classification"""

    classification_type = ClassificationType.AMPLIFICATION
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.FEATURE_AMPLIFICATION
