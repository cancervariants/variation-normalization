"""The module for Amino Acid Deletion Validation."""
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import AminoAcidDeletionToken
from typing import List, Optional
from variation.validators.validator import Validator
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.token_response_schema import Token
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.data_sources import SeqRepoAccess, TranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from .amino_acid_base import AminoAcidBase
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class AminoAcidDeletion(Validator):
    """The Amino Acid Deletion Validator class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTA,
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
            uta
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
                                  normalize_endpoint) -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of Tokens
        :param list transcripts: A list of transcript strings
        :param Classification classification: A classification for a list of
            tokens
        :param list results: A list to store validation result objects
        :param list gene_tokens: List of GeneMatchTokens
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        """
        valid_alleles = list()
        if 'HGVS' in classification.matching_tokens:
            is_hgvs = True
        else:
            is_hgvs = False

        mane_data = {
            'mane_select': dict(),
            'mane_plus_clinical': dict(),
            'longest_compatible_remaining': dict()
        }

        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                allele, t, hgvs_expr, is_ensembl = \
                    self.get_allele_with_context(classification, t, s, errors)

                if not allele:
                    errors.append("Unable to find allele.")
                else:
                    if len(allele['state']['sequence']) == 3:
                        allele['state']['sequence'] = \
                            self._amino_acid_cache.convert_three_to_one(
                                allele['state']['sequence'])

                # Check ref start/end protein matches expected
                if not errors:
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
                    if mane:
                        prefix = \
                            f"{mane['refseq']}:{s.reference_sequence.lower()}."
                        dels = f"{s.start_aa_del}{mane['pos'][0]}"
                        if s.start_pos_del is not None and s.end_pos_del is not None:  # noqa: E501
                            dels += f"_{s.end_aa_del}{mane['pos'][1]}"
                        mane_hgvs_expr = f"{prefix}{dels}del"
                        self.add_mane_data(mane_hgvs_expr, mane, mane_data, s)

                self.add_validation_result(allele, valid_alleles, results,
                                           classification, s, t, gene_tokens,
                                           errors)

                if is_hgvs:
                    break

        self.add_mane_to_validation_results(
            mane_data, valid_alleles, results, classification, gene_tokens
        )

    def get_hgvs_expr(self, classification, t, s, is_hgvs) -> str:
        """Return HGVS expression

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :param Token s: The classification token
        :param bool is_hgvs: Whether or not classification is HGVS token
        :return: HGVS expression
        """
        if not is_hgvs:
            prefix = f"{t}:{s.reference_sequence.lower()}."
            dels = f"{s.start_aa_del}{s.start_pos_del}"
            if s.start_pos_del is not None and s.end_pos_del is not None:
                dels += f"_{s.end_aa_del}{s.end_pos_del}"
            hgvs_expr = f"{prefix}{dels}del"
        else:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type == 'HGVS']
            if not hgvs_token:
                hgvs_token = \
                    [t for t in classification.all_tokens if
                     isinstance(t, Token) and t.token_type == 'ReferenceSequence']  # noqa: E501
            hgvs_token = hgvs_token[0]
            hgvs_expr = hgvs_token.input_string
        return hgvs_expr

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
