"""The module for Genomic Deletion Range Validation."""
from .validator import Validator
from variation.schemas.classification_response_schema import \
    ClassificationType, Classification
from variation.schemas.token_response_schema import \
    GenomicDeletionRangeToken, Token
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
from gene.query import QueryHandler as GeneQueryHandler
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicDeletionRange(Validator):
    """The Genomic Deletion Range Validator class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTA, dp: SeqRepoDataProxy, tlr: Translator,
                 gene_normalizer: GeneQueryHandler):
        """Initialize the Genomic Deletion Range validator.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: Gene symbol tokenizer
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param UTA uta: Access to UTA queries
        :param GeneQueryHandler gene_normalizer: Access to gene-normalizer
        """
        super().__init__(
            seq_repo_access, transcript_mappings, gene_symbol, mane_transcript,
            uta, dp, tlr, gene_normalizer
        )
        self.hgvs_dup_del_mode = HGVSDupDelMode(seq_repo_access)

    def get_transcripts(self, gene_tokens: List,
                        classification: Classification,
                        errors: List) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        return self.get_genomic_transcripts(classification, errors)

    def get_valid_invalid_results(
            self, classification_tokens: List, transcripts: List,
            classification: Classification, results: List, gene_tokens: List,
            normalize_endpoint: bool, mane_data_found: Dict,
            is_identifier: bool, hgvs_dup_del_mode: HGVSDupDelModeEnum
    ) -> None:
        """Add validation result objects to a list of results.

        :param List classification_tokens: A list of classification Tokens
        :param List transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param List results: Stores validation result objects
        :param List gene_tokens: List of GeneMatchTokens for a classification
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param Dict mane_data_found: MANE Transcript information found
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

                result = self._get_variation(s, t, errors, gene_tokens,
                                             hgvs_dup_del_mode)
                variation = result['variation']

                if not errors:
                    self._get_normalize_variation(
                        gene_tokens, s, t, errors, hgvs_dup_del_mode,
                        mane_data_found)

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

    def _get_variation(
            self, s: Token, t: str, errors: List, gene_tokens: List,
            hgvs_dup_del_mode: HGVSDupDelModeEnum) -> Optional[Dict]:
        """Get variation data.

        :param Token s: Classification token
        :param str t: Accession
        :param List errors: List of errors
        :param List gene_tokens: List of GeneMatchTokens for a classification
        :param HGVSDupDelMode hgvs_dup_del_mode: Mode to use for interpreting
            HGVS duplications and deletions
        :return: Dictionary containing start/end position changes and variation
        """
        variation, start, end = None, None, None
        ival, grch38 = self._get_ival(t, s, errors, gene_tokens)

        if not errors:
            if grch38:
                t = grch38['ac']

            allele = self.to_vrs_allele_ranges(
                t, s.reference_sequence, s.alt_type, errors, ival)

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

    def _get_normalize_variation(
            self, gene_tokens: List, s: Token, t: str, errors: List,
            hgvs_dup_del_mode: HGVSDupDelModeEnum,
            mane_data_found: Dict) -> None:
        """Get variation that will be returned in normalize endpoint.

        :param List gene_tokens: List of gene tokens
        :param Token s: Classification token
        :param str t: Accession
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Mode to use for
            interpreting HGVS duplications and deletions
        :param Dict mane_data_found: MANE Transcript data found for given query
        """
        ival, grch38 = self._get_ival(t, s, errors, gene_tokens, is_norm=True)
        if not errors:
            if grch38:
                t = grch38['ac']

            allele = self.to_vrs_allele_ranges(
                t, s.reference_sequence, s.alt_type, errors, ival)

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

    def _get_ival(
            self, t: str, s: Token, errors: List, gene_tokens: List,
            is_norm: bool = False
    ) -> Optional[Tuple[Optional[models.SequenceInterval], Optional[Dict]]]:
        """Get ival for variations with ranges.

        :param str t: Accession
        :param Token t: Classification token
        :param List errors: List of errors
        :param List gene_tokens: LIst of gene tokens
        :param bool is_norm: `True` if normalize endpoint is being used.
            `False` otherwise.
        :return: Sequence Interval and GRCh38 data if normalize endpoint
            is being used
        """
        ival, grch38, start1, start2, end1, end2 = None, None, None, None, None, None  # noqa: E501
        # (#_#)_(#_#)
        if is_norm:
            grch38_start = self.mane_transcript.g_to_grch38(
                t, s.start_pos1_del, s.start_pos2_del
            )
            if grch38_start:
                start1, start2 = grch38_start['pos']

                grch38_end = self.mane_transcript.g_to_grch38(
                    t, s.end_pos1_del, s.end_pos2_del
                )
                t = grch38_start['ac']
                if grch38_end:
                    end1, end2 = grch38_end['pos']
                    grch38 = grch38_end  # Pos doesn't really matter in return
                else:
                    return ival, grch38
            else:
                return ival, grch38
        else:
            start1 = s.start_pos1_del
            start2 = s.start_pos2_del
            end1 = s.end_pos1_del
            end2 = s.end_pos2_del

        if gene_tokens:
            self._validate_gene_pos(
                gene_tokens[0].token, t, start1, start2, errors,
                pos3=end1, pos4=end2
            )
        else:
            for pos in [start1, start2, end1, end2]:
                self._check_index(t, pos, errors)

        if not errors and start1 and start2 and end1 and end2:
            ival = self._get_ival_certain_range(
                start1, start2, end1, end2
            )

        return ival, grch38

    def get_gene_tokens(
            self, classification: Classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self) -> str:
        """Return the variation name."""
        return 'genomic deletion range'

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Genomic Deletion Range.

        :param Token t: Classification token
        """
        return t.token_type == 'GenomicDeletionRange'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Genomic Deletion Range.

        :param ClassificationType classification_type: Classification type
        :return: `True` if classification type matches, `False` otherwise
        """
        return classification_type == \
            ClassificationType.GENOMIC_DELETION_RANGE

    def human_description(self, transcript,
                          token: GenomicDeletionRangeToken) -> str:
        """Return a human description of the identified variation."""
        descr = f"A Genomic Deletion from" \
                f" ({token.start_pos1_del}_{token.start_pos2_del}) to " \
                f"({token.end_pos1_del}_{token.end_pos2_del}) on {transcript}"
        return descr

    def concise_description(self, transcript, token):
        """Return a concise description of the identified variation."""
        return f"{transcript}:g.({token.start_pos1_del}_" \
               f"{token.start_pos2_del})_({token.end_pos1_del}_" \
               f"{token.end_pos2_del})"
