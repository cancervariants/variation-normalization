"""The module for Amino Acid Insertion Validation."""
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import AminoAcidInsertionToken
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


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class AminoAcidInsertion(Validator):
    """The Amino Acid Insertion Validator class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTA, dp: SeqRepoDataProxy, tlr: Translator,
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
        :param amino_acid_cache: Amino Acid codes and conversions
        """
        super().__init__(
            seq_repo_access, transcript_mappings, gene_symbol, mane_transcript,
            uta, dp, tlr
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
                                  is_identifier) -> None:
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
        """
        valid_alleles = list()
        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                t = self.get_accession(t, classification)
                allele = self.to_vrs_allele(
                    t, s.start_pos_flank, s.end_pos_flank,
                    s.reference_sequence, s.alt_type, errors,
                    alt=s.inserted_sequence
                )

                if not errors:
                    # Check ref start/end protein matches expected
                    self.amino_acid_base.check_ref_aa(
                        t, s.start_aa_flank, s.start_pos_flank, errors
                    )

                    if not errors:
                        self.amino_acid_base.check_ref_aa(
                            t, s.end_aa_flank, s.end_pos_flank, errors
                        )

                if not errors:
                    mane = self.mane_transcript.get_mane_transcript(
                        t, s.start_pos_flank, s.end_pos_flank,
                        s.reference_sequence,
                        normalize_endpoint=normalize_endpoint
                    )
                    self.add_mane_data(mane, mane_data_found,
                                       s.reference_sequence, s.alt_type, s,
                                       gene_tokens, alt=s.inserted_sequence)

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
        return 'amino acid insertion'

    def is_token_instance(self, t):
        """Check that token is Amino Acid Insertion."""
        return t.token_type == 'AminoAcidInsertion'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Amino Acid Insertion.
        """
        return classification_type == ClassificationType.AMINO_ACID_INSERTION

    def concise_description(self, transcript, token) -> str:
        """Return a HGVS description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression
        """
        return f'{transcript}:p.' \
               f'{token.start_aa_flank}{token.start_pos_flank}_' \
               f'{token.end_aa_flank}{token.end_pos_flank}' \
               f'ins{token.inserted_sequence}'

    def human_description(self, transcript,
                          token: AminoAcidInsertionToken) -> str:
        """Return a human description of the identified variation."""
        return f"The insertion of amino acid(s) {token.inserted_sequence} " \
               f"between amino acids " \
               f"{token.start_aa_flank}{token.start_pos_flank} and " \
               f"{token.end_aa_flank}{token.end_pos_flank}"
