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
from typing import List


class Validate:
    """The validation class."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTA,
                 amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the validate class.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: Gene symbol tokenizer
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param UTA uta: Access to UTA queries
        :param amino_acid_cache: Amino Acid codes and conversions
        """
        self.validators = [
            AminoAcidSubstitution(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta, amino_acid_cache
            ),
            PolypeptideTruncation(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta, amino_acid_cache
            ),
            SilentMutation(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta, amino_acid_cache
            ),
            CodingDNASubstitution(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta
            ),
            GenomicSubstitution(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta
            ),
            CodingDNASilentMutation(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta
            ),
            GenomicSilentMutation(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta
            ),
            AminoAcidDelIns(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta, amino_acid_cache
            ),
            CodingDNADelIns(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta
            ),
            GenomicDelIns(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta
            ),
            AminoAcidDeletion(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta, amino_acid_cache
            ),
            CodingDNADeletion(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta
            ),
            GenomicDeletion(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta
            ),
            AminoAcidInsertion(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta, amino_acid_cache
            ),
            CodingDNAInsertion(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta
            ),
            GenomicInsertion(
                seqrepo_access, transcript_mappings, gene_symbol,
                mane_transcript, uta
            )
        ]

    def perform(self, classifications: List[Classification],
                normalize_endpoint, warnings=None) \
            -> ValidationSummary:
        """Validate a list of classifications."""
        valid_possibilities = list()
        invalid_possibilities = list()
        if not warnings:
            warnings = list()

        found_classification = False
        for classification in classifications:
            for validator in self.validators:
                if validator.validates_classification_type(
                        classification.classification_type):
                    results = validator.validate(classification,
                                                 normalize_endpoint)
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
