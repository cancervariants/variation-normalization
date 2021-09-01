"""The module for Genomic Uncertain Deletion Validation."""
from .validator import Validator
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import \
    GenomicUncertainDeletionToken
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicUncertainDeletion(Validator):
    """The Genomic UncertainDeletion Validator class."""

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
                                  is_identifier) -> None:
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
        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                t = self.get_accession(t, classification)

                allele = self.to_vrs_allele(t, s.start_pos2_del,
                                            s.end_pos1_del,
                                            s.reference_sequence, s.alt_type,
                                            errors)
                cnv = self.to_vrs_cnv(t, allele, 'del')
                if not cnv:
                    errors.append(f"Unable to get CNV for {t}")

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
                        s.alt_type, s, gene_tokens
                    )

                self.add_validation_result(
                    cnv, valid_alleles, results,
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
