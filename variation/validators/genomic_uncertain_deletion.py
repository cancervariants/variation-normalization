"""The module for Genomic Uncertain Deletion Validation."""
from .validator import Validator
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import \
    GenomicUncertainDeletionToken
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
import logging
from variation.schemas.normalize_response_schema \
    import HGVSDupDelMode as HGVSDupDelModeEnum
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

                ival = models.SequenceInterval(
                    start=models.IndefiniteRange(
                        value=s.start_pos2_del - 1,
                        comparator="<="
                    ),
                    end=models.IndefiniteRange(
                        value=s.end_pos1_del,
                        comparator=">="
                    )
                )
                allele = self.to_vrs_allele_ranges(
                    s, t, s.reference_sequence, s.alt_type, errors,
                    ival
                )
                if hgvs_dup_del_mode == HGVSDupDelModeEnum.DEFAULT or \
                        hgvs_dup_del_mode == HGVSDupDelMode.CNV:
                    variation = self.to_vrs_cnv(t, allele, 'del',
                                                uncertain=True)
                    if not variation:
                        errors.append(f"Unable to get CNV for {t}")
                elif hgvs_dup_del_mode == HGVSDupDelModeEnum.LITERAL_SEQ_EXPR:
                    variation = allele
                    if not variation:
                        errors.append(f"Unable to get allele for {t}")
                else:
                    variation = None
                    errors.append("Genomic uncertain deletions cannot have"
                                  " repeated seq expressions")

                if not errors:
                    grch38 = self.mane_transcript.g_to_grch38(
                        t, s.start_pos2_del, s.end_pos1_del)

                    if grch38:
                        mane = dict(
                            gene=None,
                            refseq=grch38['ac'] if grch38['ac'].startswith('NC') else None,  # noqa: E501
                            ensembl=grch38['ac'] if grch38['ac'].startswith('ENSG') else None,  # noqa: E501
                            pos=grch38['pos'],
                            strand=None,
                            status='GRCh38'
                        )
                    else:
                        mane = None

                    self.add_mane_data(
                        mane, mane_data_found, s.reference_sequence,
                        s.alt_type, s, gene_tokens,
                        hgvs_dup_del_mode=hgvs_dup_del_mode
                    )

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
