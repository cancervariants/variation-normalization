"""The module for Coding DNA Substitution Validation."""
from .single_nucleotide_variation_base import SingleNucleotideVariationBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import CodingDNASubstitutionToken
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
import logging

# TODO:
#  LRG_ (LRG_199t1:c)


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class CodingDNASubstitution(SingleNucleotideVariationBase):
    """The Coding DNA Substitution Validator class."""

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

                    allele = self.to_vrs_allele(t, s.position, s.position,
                                                s.reference_sequence,
                                                s.alt_type, errors,
                                                cds_start=cds_start,
                                                alt=s.new_nucleotide)
                else:
                    allele = None
                    cds_start = 0
                    errors.append(f"Unable to get CDS start for {t}")

                if not errors:
                    ref_nuc = self.seqrepo_access.get_sequence(
                        t, s.position + cds_start
                    )
                    self.check_ref_nucleotide(ref_nuc, s.ref_nucleotide,
                                              s.position, t, errors)

                if not errors:
                    mane = self.mane_transcript.get_mane_transcript(
                        t, s.position, s.position, s.reference_sequence,
                        ref=s.ref_nucleotide,
                        normalize_endpoint=normalize_endpoint
                    )

                    self.add_mane_data(
                        mane, mane_data_found, s.reference_sequence,
                        s.alt_type, s, gene_tokens,
                        alt=s.new_nucleotide
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
        return 'coding dna substitution'

    def is_token_instance(self, t):
        """Check that token is Coding DNA Substitution."""
        return t.token_type == 'CodingDNASubstitution'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is coding dna
        substitution.
        """
        return classification_type == ClassificationType.CODING_DNA_SUBSTITUTION  # noqa: E501

    def human_description(self, transcript,
                          token: CodingDNASubstitutionToken) -> str:
        """Return a human description of the identified variation."""
        return f'A coding DNA substitution from {token.ref_nucleotide}' \
               f' to {token.new_nucleotide} at position ' \
               f'{token.position} on transcript {transcript}'
