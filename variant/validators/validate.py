"""Module for Validation."""
from variant.schemas.validation_response_schema import ValidationSummary
from variant.schemas.classification_response_schema import Classification
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import AminoAcidCache
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
from typing import List


class Validate:
    """The validation class."""

    def __init__(self, seq_repo_client: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the validate class."""
        self.validators = [
            AminoAcidSubstitution(seq_repo_client, transcript_mappings,
                                  gene_symbol, amino_acid_cache),
            PolypeptideTruncation(seq_repo_client, transcript_mappings,
                                  gene_symbol, amino_acid_cache),
            SilentMutation(seq_repo_client, transcript_mappings,
                           gene_symbol, amino_acid_cache),
            CodingDNASubstitution(seq_repo_client, transcript_mappings,
                                  gene_symbol),
            GenomicSubstitution(seq_repo_client, transcript_mappings,
                                gene_symbol),
            CodingDNASilentMutation(seq_repo_client, transcript_mappings,
                                    gene_symbol),
            GenomicSilentMutation(seq_repo_client, transcript_mappings,
                                  gene_symbol),
            AminoAcidDelIns(seq_repo_client, transcript_mappings, gene_symbol,
                            amino_acid_cache),
            CodingDNADelIns(seq_repo_client, transcript_mappings, gene_symbol),
            GenomicDelIns(seq_repo_client, transcript_mappings, gene_symbol),
            AminoAcidDeletion(seq_repo_client, transcript_mappings,
                              gene_symbol, amino_acid_cache),
            CodingDNADeletion(seq_repo_client, transcript_mappings,
                              gene_symbol),
            GenomicDeletion(seq_repo_client, transcript_mappings, gene_symbol)
        ]

    def perform(self, classifications: List[Classification]) \
            -> ValidationSummary:
        """Validate a list of classifications."""
        valid_possibilities = list()
        invalid_possibilities = list()

        for classification in classifications:
            for validator in self.validators:
                if validator.validates_classification_type(
                        classification.classification_type):
                    results = validator.validate(classification)
                    for res in results:
                        if res.is_valid:
                            valid_possibilities.append(res)
                        else:
                            invalid_possibilities.append(res)

        return ValidationSummary(
            valid_results=valid_possibilities,
            invalid_results=invalid_possibilities
        )
