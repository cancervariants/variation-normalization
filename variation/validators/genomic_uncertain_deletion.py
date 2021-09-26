"""The module for Genomic Uncertain Deletion Validation."""
from .validator import Validator
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import \
    GenomicUncertainDeletionToken
from typing import List, Optional, Dict, Tuple
from variation.schemas.token_response_schema import GeneMatchToken
import logging
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.data_sources import SeqRepoAccess, TranscriptMappings, UTA
from variation.tokenizers import GeneSymbol
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from ga4gh.vrs import models

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicUncertainDeletion(Validator):
    """The Genomic UncertainDeletion Validator class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTA, dp: SeqRepoDataProxy, tlr: Translator):
        """Initialize the Genomic Uncertain Deletion validator.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: Gene symbol tokenizer
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param UTA uta: Access to UTA queries
        """
        super().__init__(
            seq_repo_access, transcript_mappings, gene_symbol, mane_transcript,
            uta, dp, tlr
        )
        self.hgvs_dup_del_mode = HGVSDupDelMode(seq_repo_access)

    def get_transcripts(self, gene_tokens, classification, errors)\
            -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param list gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param list errors: List of errors
        :return: List of transcript accessions
        """
        return self.get_genomic_transcripts(classification, errors)

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
        :param str hgvs_dup_del_mode: Must be: `default`, `cnv`,
            `repeated_seq_expr`, `literal_seq_expr`.
            This parameter determines how to represent HGVS dup/del expressions
            as VRS objects.
        """
        valid_alleles = list()
        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                t = self.get_accession(t, classification)

                result = self._get_variation(s, t, errors, hgvs_dup_del_mode)
                variation = result['variation']
                start = result['start']
                end = result['end']

                if not errors:
                    self._get_normalize_variation(
                        gene_tokens, s, t, errors, hgvs_dup_del_mode,
                        mane_data_found, start, end)

                self.add_validation_result(
                    variation, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

                if is_identifier:
                    break

        self.add_mane_to_validation_results(
            mane_data_found, valid_alleles, results,
            classification, gene_tokens
        )

    def _get_variation(self, s, t, errors, hgvs_dup_del_mode)\
            -> Optional[Dict]:
        """Get variation data.

        :param Token s: Classification token
        :param str t: Accession
        :param list errors: List of errors
        :param HGVSDupDelMode hgvs_dup_del_mode: Mode to use for interpreting
            HGVS duplications and deletions
        :return: Dictionary containing start/end position changes and variation
        """
        variation, start, end = None, None, None
        ival, _ = self._get_ival(t, s, errors)

        allele = self.to_vrs_allele_ranges(
            s, t, s.reference_sequence, s.alt_type,
            errors, ival
        )
        if start is not None and end is not None:
            pos = (start, end)
        else:
            pos = None
        variation = self.hgvs_dup_del_mode.interpret_variation(
            t, s.alt_type, allele, errors,
            hgvs_dup_del_mode, pos=pos)

        return {
            'start': start,
            'end': end,
            'variation': variation
        }

    def _get_normalize_variation(self, gene_tokens, s, t, errors,
                                 hgvs_dup_del_mode, mane_data_found, start,
                                 end) -> None:
        """Get variation that will be returned in normalize endpoint.

        :param list gene_tokens: List of gene tokens
        :param Token s: Classification token
        :param str t: Accession
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Mode to use for
            interpreting HGVS duplications and deletions
        :param dict mane_data_found: MANE Transcript data found for given query
        :param int start: Start pos change
        :param int end: End pos change
        """
        if not gene_tokens:
            ival, grch38 = self._get_ival(t, s, errors, is_norm=True)

            if not errors:
                allele = self.to_vrs_allele_ranges(
                    s, t, s.reference_sequence,
                    s.alt_type, errors, ival)

                grch38_variation = \
                    self.hgvs_dup_del_mode.interpret_variation(
                        t, s.alt_type, allele, errors,
                        hgvs_dup_del_mode
                    )

                if grch38_variation:
                    self._add_dict_to_mane_data(
                        grch38['ac'], s, grch38_variation,
                        mane_data_found, 'GRCh38'
                    )

    def _get_ival(self, t, s, errors, is_norm=False)\
            -> Optional[Tuple[models.SequenceInterval, Dict]]:
        """Get ival for variations with ranges.

        :param str t: Accession
        :param Token t: Classification token
        :param list errors: List of errors
        :param bool is_norm: `True` if normalize endpoint is being used.
            `False` otherwise.
        :return: Sequence Interval and GRCh38 data if normalize endpoint
            is being used
        """
        ival = None
        start = None
        end = None
        grch38 = None
        if s.start_pos1_del == '?' and s.end_pos2_del == '?':
            # format: (?_#)_(#_?)
            if is_norm:
                grch38 = self.mane_transcript.g_to_grch38(
                    t, s.start_pos2_del, s.end_pos1_del
                )
                if grch38:
                    start, end = grch38['pos']
            else:
                start = s.start_pos2_del
                end = s.end_pos1_del

            if start and end:
                ival = models.SequenceInterval(
                    start=self._get_start_indef_range(start),
                    end=self._get_end_indef_range(end)
                )
        return ival, grch38

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
        return 'genomic uncertain deletion'

    def is_token_instance(self, t):
        """Check that token is Genomic Uncertain Deletion."""
        return t.token_type == 'GenomicUncertainDeletion'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Genomic Uncertain Deletion.
        """
        return classification_type == \
            ClassificationType.GENOMIC_UNCERTAIN_DELETION

    def human_description(self, transcript,
                          token: GenomicUncertainDeletionToken) -> str:
        """Return a human description of the identified variation."""
        descr = f"A Genomic Uncertain Deletion from" \
                f" (?_{token.start_pos2_del}) to {token.end_pos1_del}_? " \
                f"on {transcript}"
        return descr

    def concise_description(self, transcript, token):
        """Return a consice description of the identified variation."""
        return f"{transcript}:g.(?_{token.start_pos2_del})_" \
               f"({token.end_pos1_del}_?)"
