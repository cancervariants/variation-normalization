"""The module for Coding DNA DelIns Validation."""
from variation.validators.delins_base import DelInsBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import CodingDNADelInsToken
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class CodingDNADelIns(DelInsBase):
    """The Coding DNA DelIns Validator class."""

    def get_transcripts(self, gene_tokens, classification, errors)\
            -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param list gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param list errors: List of errors
        :return: List of transcript accessions
        """
        return self.get_coding_dna_transcripts(gene_tokens, errors)

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

                cds_start_end = self.uta.get_cds_start_end(t)
                if cds_start_end is not None:
                    cds_start = cds_start_end[0]
                    allele = self.to_vrs_allele(
                        t, s.start_pos_del, s.end_pos_del,
                        s.reference_sequence, s.alt_type,
                        errors, cds_start=cds_start, alt=s.inserted_sequence1
                    )
                else:
                    allele = None
                    errors.append(f"Unable to find CDS start for {t}")

                if not errors:
                    self.check_pos_index(t, s, errors)

                if not errors:
                    mane = self.mane_transcript.get_mane_transcript(
                        t, s.start_pos_del, s.end_pos_del,
                        s.reference_sequence,
                        normalize_endpoint=normalize_endpoint
                    )

                    self.add_mane_data(
                        mane, mane_data_found, s.reference_sequence,
                        s.alt_type, s, gene_tokens, alt=s.inserted_sequence1
                    )

                self.add_validation_result(
                    allele, valid_alleles, results,
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
        return self.get_coding_dna_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
        return 'coding dna delins'

    def is_token_instance(self, t):
        """Check that token is Coding DNA DelIns."""
        return t.token_type == 'CodingDNADelIns'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Coding DNA DelIns.
        """
        return classification_type == ClassificationType.CODING_DNA_DELINS

    def human_description(self, transcript,
                          token: CodingDNADelInsToken) -> str:
        """Return a human description of the identified variation."""
        if token.start_pos_del is not None and token.end_pos_del is not None:
            position = f"{token.start_pos_del} to {token.end_pos_del}"
        else:
            position = token.start_pos_del

        if token.inserted_sequence1 is not None and \
                token.inserted_sequence2 is not None:
            sequence = f"{token.inserted_sequence1} to " \
                       f"{token.inserted_sequence2}"
        else:
            sequence = token.inserted_sequence1

        return f"A Coding DNA DelIns deletion of {position} replaced by " \
               f"{sequence} on transcript {transcript}"
