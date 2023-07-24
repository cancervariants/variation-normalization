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

    PROTEIN_SUBSTITUTION = "protein substitution"
    PROTEIN_STOP_GAIN = "protein stop gain"
    PROTEIN_REFERENCE_AGREE = "reference agree"
    PROTEIN_DELINS = "protein delins"
    CDNA_SUBSTITUTION = "cdna substitution"
    GENOMIC_SUBSTITUTION = "genomic substitution"
    CDNA_REFERENCE_AGREE = "cdna reference agree"
    GENOMIC_REFERENCE_AGREE = "genomic reference agree"
    CDNA_DELINS = "cdna delins"
    GENOMIC_DELINS = "genomic delins"
    PROTEIN_DELETION = "protein deletion"
    CDNA_DELETION = "cdna deletion"
    GENOMIC_DELETION = "genomic deletion"
    GENOMIC_DELETION_AMBIGUOUS = "genomic deletion ambiguous"
    PROTEIN_INSERTION = "protein insertion"
    CDNA_INSERTION = "cdna insertion"
    GENOMIC_INSERTION = "genomic insertion"
    GENOMIC_DUPLICATION = "genomic duplication"
    GENOMIC_DUPLICATION_AMBIGUOUS = "genomic duplication ambiguous"
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

    classification_type = ClassificationType.PROTEIN_SUBSTITUTION
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.PROTEIN_SUBSTITUTION


class GenomicSubstitutionClassification(Classification, Substitution):

    classification_type = ClassificationType.GENOMIC_SUBSTITUTION
    molecule_context = MoleculeContext.GENOMIC
    so_id: Union[Literal[SequenceOntology.SNV], Literal[SequenceOntology.MNV]]


class CdnaSubstitutionClassification(Classification, Substitution):

    classification_type = ClassificationType.CDNA_SUBSTITUTION
    molecule_context = MoleculeContext.TRANSCRIPT
    so_id: Union[Literal[SequenceOntology.SNV], Literal[SequenceOntology.MNV]]


class ProteinStopGainClassification(Classification, StopGain):

    classification_type = ClassificationType.PROTEIN_STOP_GAIN
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.STOP_GAIN


class ProteinReferenceAgreeClassification(Classification, ProteinReferenceAgree):

    classification_type = ClassificationType.PROTEIN_REFERENCE_AGREE
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.NO_SEQUENCE_ALTERATION


class CdnaReferenceAgreeClassification(Classification, ReferenceAgree):

    classification_type = ClassificationType.CDNA_REFERENCE_AGREE
    molecule_context = MoleculeContext.TRANSCRIPT
    so_id = SequenceOntology.NO_SEQUENCE_ALTERATION


class GenomicReferenceAgreeClassification(Classification, ReferenceAgree):

    classification_type = ClassificationType.GENOMIC_REFERENCE_AGREE
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.NO_SEQUENCE_ALTERATION


class ProteinInsertionClassification(Classification, ProteinInsertion):

    classification_type = ClassificationType.PROTEIN_INSERTION
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.PROTEIN_INSERTION


class CdnaInsertionClassification(Classification, Insertion):

    classification_type = ClassificationType.CDNA_INSERTION
    molecule_context = MoleculeContext.TRANSCRIPT
    so_id = SequenceOntology.INSERTION


class GenomicInsertionClassification(Classification, Insertion):

    classification_type = ClassificationType.GENOMIC_INSERTION
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.INSERTION


class ProteinDeletionClassification(Classification, ProteinDeletion):

    classification_type = ClassificationType.PROTEIN_DELETION
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.PROTEIN_DELETION


class GenomicDeletionClassification(Classification, Deletion):

    classification_type = ClassificationType.GENOMIC_DELETION
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.DELETION


class CdnaDeletionClassification(Classification, Deletion):

    classification_type = ClassificationType.CDNA_DELETION
    molecule_context = MoleculeContext.TRANSCRIPT
    so_id = SequenceOntology.DELETION


class ProteinDelInsClassification(Classification, ProteinDelIns):

    classification_type = ClassificationType.PROTEIN_DELINS
    molecule_context = MoleculeContext.PROTEIN
    so_id = SequenceOntology.DELINS


class CdnaDelInsClassification(Classification, DelIns):

    classification_type = ClassificationType.CDNA_DELINS
    molecule_context = MoleculeContext.TRANSCRIPT
    so_id = SequenceOntology.DELINS


class GenomicDelInsClassification(Classification, DelIns):

    classification_type = ClassificationType.GENOMIC_DELINS
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.DELINS


class GenomicDuplicationClassification(Classification, Duplication):

    classification_type = ClassificationType.GENOMIC_DUPLICATION
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.DUPLICATION


class AmbiguousType(str, Enum):
    """Helps determine the kind of ambiguous variant"""

    AMBIGUOUS_1 = "(#_#)_(#_#)"
    AMBIGUOUS_2 = "(?_#)_(#_?)"
    AMBIGUOUS_3 = "(#_?)_(?_#)"  # Not yet supported
    AMBIGUOUS_4 = "(#_#)_#"  # Not yet supported
    AMBIGUOUS_5 = "(?_#)_#"
    AMBIGUOUS_6 = "#_(#_#)"  # Not yet supported
    AMBIGUOUS_7 = "#_(#_?)"


class GenomicDuplicationAmbiguousClassification(Classification, DupDelAmbiguous):

    classification_type = ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.DUPLICATION
    ambiguous_type: AmbiguousType


class GenomicDeletionAmbiguousClassification(Classification, DupDelAmbiguous):

    classification_type = ClassificationType.GENOMIC_DELETION_AMBIGUOUS
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.DELETION
    ambiguous_type: AmbiguousType


class AmplificationClassification(Classification):

    classification_type = ClassificationType.AMPLIFICATION
    molecule_context = MoleculeContext.GENOMIC
    so_id = SequenceOntology.FEATURE_AMPLIFICATION
