"""The module for Coding DNA Substitution Validation."""
from .single_nucleotide_variation_base import SingleNucleotideVariationBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import CodingDNASubstitutionToken
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.token_response_schema import Token
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

    def get_hgvs_expr(self, classification, t, s, is_hgvs) -> str:
        """Return HGVS expression

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :param Token s: The classification token
         :param bool is_hgvs: Whether or not classification is HGVS token
        :return: hgvs expression
        """
        hgvs_from_transcript = f"{t}:{s.reference_sequence.lower()}." \
                               f"{s.position}{s.ref_nucleotide}" \
                               f">{s.new_nucleotide}"
        if not is_hgvs:
            hgvs_expr = hgvs_from_transcript
        else:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type
                          in ['HGVS', 'ReferenceSequence']][0]
            hgvs_expr = hgvs_token.input_string
        return hgvs_expr

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

                mane = self.mane_transcript.get_mane_transcript(
                    t, s.position, s.position, s.reference_sequence,
                    ref=s.ref_nucleotide, normalize_endpoint=normalize_endpoint
                )

                if mane:
                    if 'coding_start_site' in mane.keys():
                        ref = self.seqrepo_access.sequence_at_position(
                            mane['refseq'],
                            mane['pos'][0] + mane['coding_start_site']
                        )

                        mane_hgvs_expr = f"{mane['refseq']}:" \
                                         f"{s.reference_sequence.lower()}." \
                                         f"{mane['pos'][0]}{ref}" \
                                         f">{s.new_nucleotide}"
                        self.add_mane_data(mane_hgvs_expr, mane, mane_data, s)
                    else:
                        errors.append("No coding start site found.")

                if not allele:
                    errors.append("Unable to find allele.")
                else:
                    ref_nuc = \
                        self.seqrepo_access.sequence_at_position(t, s.position)
                    self.check_ref_nucleotide(ref_nuc, s, t, errors)

                self.add_validation_result(
                    allele, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

                if is_hgvs:
                    break

        self.add_mane_to_validation_results(
            mane_data, valid_alleles, results, classification, gene_tokens
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
