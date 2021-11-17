"""Module for Validation."""
from variation.schemas.validation_response_schema import ValidationSummary
from variation.schemas.classification_response_schema import Classification
from variation.data_sources import TranscriptMappings, SeqRepoAccess, UTA
from variation.mane_transcript import MANETranscript
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from .amino_acid_substitution import AminoAcidSubstitution
from .polypeptide_truncation import PolypeptideTruncation
from .silent_mutation import SilentMutation
from .coding_dna_substitution import CodingDNASubstitution
from .coding_dna_silent_mutation import CodingDNASilentMutation
from .genomic_silent_mutation import GenomicSilentMutation
from .genomic_substitution import GenomicSubstitution
from .amino_acid_delins import AminoAcidDelIns
from .coding_dna_delins import CodingDNADelIns
from .genomic_delins import GenomicDelIns
from .amino_acid_deletion import AminoAcidDeletion
from .coding_dna_deletion import CodingDNADeletion
from .genomic_deletion import GenomicDeletion
from .amino_acid_insertion import AminoAcidInsertion
from .coding_dna_insertion import CodingDNAInsertion
from .genomic_insertion import GenomicInsertion
from .genomic_uncertain_deletion import GenomicUncertainDeletion
from .genomic_duplication import GenomicDuplication
from .genomic_deletion_range import GenomicDeletionRange
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from typing import List
from gene.query import QueryHandler as GeneQueryHandler
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


class Validate:
    """The validation class."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTA, dp: SeqRepoDataProxy, tlr: Translator,
                 amino_acid_cache: AminoAcidCache,
                 gene_normalizer: GeneQueryHandler) -> None:
        """Initialize the validate class.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: Gene symbol tokenizer
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param UTA uta: Access to UTA queries
        :param Translator tlr: Translator class
        :param GeneQueryHandler gene_normalizer: Access to gene-normalizer
        :param amino_acid_cache: Amino Acid codes and conversions
        """
        params = [
            seqrepo_access, transcript_mappings, gene_symbol,
            mane_transcript, uta, dp, tlr, gene_normalizer
        ]
        amino_acid_params = params[:]
        amino_acid_params.append(amino_acid_cache)
        self.validators = [
            AminoAcidSubstitution(*amino_acid_params),
            PolypeptideTruncation(*amino_acid_params),
            SilentMutation(*amino_acid_params),
            CodingDNASubstitution(*params),
            GenomicSubstitution(*params),
            CodingDNASilentMutation(*params),
            GenomicSilentMutation(*params),
            AminoAcidDelIns(*amino_acid_params),
            CodingDNADelIns(*params),
            GenomicDelIns(*params),
            AminoAcidDeletion(*amino_acid_params),
            CodingDNADeletion(*params),
            GenomicDeletion(*params),
            AminoAcidInsertion(*amino_acid_params),
            CodingDNAInsertion(*params),
            GenomicInsertion(*params),
            GenomicDeletionRange(*params),
            GenomicUncertainDeletion(*params),
            GenomicDuplication(*params)
        ]

    def perform(
            self, classifications: List[Classification],
            normalize_endpoint: bool, warnings: List = None,
            hgvs_dup_del_mode: HGVSDupDelModeEnum = HGVSDupDelModeEnum.DEFAULT
    ) -> ValidationSummary:
        """Validate a list of classifications.

        :param List classifications: List of classifications
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param List warnings: List of warnings
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `cnv`,
            `repeated_seq_expr`, `literal_seq_expr`.
            This parameter determines how to represent HGVS dup/del expressions
            as VRS objects.
        :return: ValidationSummary containing valid and invalid results
        """
        valid_possibilities = list()
        invalid_possibilities = list()
        if not warnings:
            warnings = list()

        found_classification = False
        for classification in classifications:
            for validator in self.validators:
                if validator.validates_classification_type(
                        classification.classification_type):
                    results = validator.validate(
                        classification, normalize_endpoint, hgvs_dup_del_mode)
                    for res in results:
                        if res.is_valid:
                            found_classification = True
                            valid_possibilities.append(res)
                        else:
                            invalid_possibilities.append(res)

                if found_classification:
                    break

        if not found_classification and not warnings:
            warnings.append("Unable find a classification")

        return ValidationSummary(
            valid_results=valid_possibilities,
            invalid_results=invalid_possibilities,
            warnings=warnings
        )
