"""The module for Amino Acid Insertion Validation."""
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import AminoAcidInsertionToken
from typing import List, Optional
from variation.validators.validator import Validator
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.token_response_schema import Token
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.data_sources import SeqRepoAccess, TranscriptMappings
from variation.mane_transcript import MANETranscript
from .amino_acid_base import AminoAcidBase
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class AminoAcidInsertion(Validator):
    """The Amino Acid Insertion Validator class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 amino_acid_cache: AminoAcidCache) \
            -> None:
        """Initialize the validator.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        """
        super().__init__(
            seq_repo_access, transcript_mappings, gene_symbol, mane_transcript
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

                    if mane:
                        refseq_ac = mane['refseq']
                        try:
                            start_aa_flank = \
                                self._amino_acid_cache.amino_acid_code_conversion[  # noqa: E501
                                    self.seqrepo_access.sequence_at_position(
                                        refseq_ac, mane['pos'][0]
                                    )
                                ]
                            end_aa_flank = \
                                self._amino_acid_cache.amino_acid_code_conversion[  # noqa: E501
                                    self.seqrepo_access.sequence_at_position(
                                        refseq_ac, mane['pos'][1]
                                    )
                                ]
                        except KeyError:
                            logger.warning(
                                f"Unable to get aa flanks for {refseq_ac} "
                                f"positions {mane['pos']}"
                            )
                        else:
                            mane_hgvs_expr = \
                                f"{refseq_ac}:{s.reference_sequence.lower()}" \
                                f".{start_aa_flank}{mane['pos'][0]}_" \
                                f"{end_aa_flank}{mane['pos'][1]}" \
                                f"ins{s.inserted_sequence}"
                            self.add_mane_data(mane_hgvs_expr, mane, mane_data,
                                               s)

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
        :return: hgvs expression
        """
        if not is_hgvs:
            hgvs_expr = f"{t}:p.{s.start_aa_flank}{s.start_pos_flank}_" \
                        f"{s.end_aa_flank}{s.end_pos_flank}" \
                        f"ins{s.inserted_sequence}"
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
