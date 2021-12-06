"""The module for Amino Acid Deletion Validation."""
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import AminoAcidDeletionToken
from typing import List, Optional
from variation.validators.validator import Validator
from variation.schemas.token_response_schema import GeneMatchToken
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.data_sources import SeqRepoAccess, TranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from .amino_acid_base import AminoAcidBase
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
import logging
from gene.query import QueryHandler as GeneQueryHandler
from variation.vrs import VRS

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class AminoAcidDeletion(Validator):
    """The Amino Acid Deletion Validator class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTA, dp: SeqRepoDataProxy, tlr: Translator,
                 gene_normalizer: GeneQueryHandler, vrs: VRS,
                 amino_acid_cache: AminoAcidCache) \
            -> None:
        """Initialize the validator.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: Gene symbol tokenizer
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param UTA uta: Access to UTA queries
        :param GeneQueryHandler gene_normalizer: Access to gene-normalizer
        :param VRS vrs: Class for creating VRS objects
        :param amino_acid_cache: Amino Acid codes and conversions
        """
        super().__init__(
            seq_repo_access, transcript_mappings, gene_symbol, mane_transcript,
            uta, dp, tlr, gene_normalizer, vrs
        )
        self._amino_acid_cache = amino_acid_cache
        self.amino_acid_base = AminoAcidBase(seq_repo_access, amino_acid_cache)
        self.mane_transcript = mane_transcript

    def get_transcripts(self, gene_tokens, classification, errors)\
            -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param list gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param list errors: List of errors
        :return: List of transcript accessions
        """
        return self.get_protein_transcripts(gene_tokens, errors)

    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results, gene_tokens,
                                  normalize_endpoint, mane_data_found,
                                  is_identifier, hgvs_dup_del_mode)\
            -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of classification Tokens
        :param list transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param list results: Stores validation result objects
        :param list gene_tokens: List of GeneMatchTokens for a classification
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param dict mane_data_found: MANE Transcript information found
        :param bool is_identifier: `True` if identifier is given for exact
            location. `False` otherwise.
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `cnv`,
            `repeated_seq_expr`, `literal_seq_expr`.
            This parameter determines how to represent HGVS dup/del expressions
            as VRS objects.
        """
        valid_alleles = list()
        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                t = self.get_accession(t, classification)
                allele = self.vrs.to_vrs_allele(
                    t, s.start_pos_del, s.end_pos_del, s.reference_sequence,
                    s.alt_type, errors)

                if not errors:
                    # Check ref start/end protein matches expected
                    self.amino_acid_base.check_ref_aa(
                        t, s.start_aa_del, s.start_pos_del, errors
                    )

                    if not errors and s.end_aa_del:
                        self.amino_acid_base.check_ref_aa(
                            t, s.end_aa_del, s.end_pos_del, errors
                        )

                if not errors:
                    mane = self.mane_transcript.get_mane_transcript(
                        t, s.start_pos_del, s.end_pos_del,
                        s.reference_sequence,
                        normalize_endpoint=normalize_endpoint
                    )
                    self.add_mane_data(
                        mane, mane_data_found, s.reference_sequence,
                        s.alt_type, s,
                    )

                self.add_validation_result(allele, valid_alleles, results,
                                           classification, s, t, gene_tokens,
                                           errors)

                if is_identifier:
                    break

        self.add_mane_to_validation_results(
            mane_data_found, valid_alleles, results,
            classification, gene_tokens
        )

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_protein_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
        return 'amino acid deletion'

    def is_token_instance(self, t):
        """Check that token is Amino Acid DelIns."""
        return t.token_type == 'AminoAcidDeletion'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Amino Acid Deletion.
        """
        return classification_type == ClassificationType.AMINO_ACID_DELETION

    def concise_description(self, transcript, token) -> str:
        """Return a HGVS description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression
        """
        dels = f"{token.start_aa_del}{token.start_pos_del}"
        if token.start_pos_del is not None and token.end_pos_del is not None:
            dels += f"_{token.end_aa_del}{token.end_pos_del}"

        return f'{transcript}:p.{dels}del'

    def human_description(self, transcript,
                          token: AminoAcidDeletionToken) -> str:
        """Return a human description of the identified variation."""
        position = f"{token.start_aa_del}{token.start_pos_del}"
        if token.start_pos_del is not None and token.end_pos_del is not None:
            position += f"to {token.end_aa_del}{token.end_pos_del}"

        return f"An Amino Acid Deletion deletion of {position} on " \
               f"transcript {transcript}"
