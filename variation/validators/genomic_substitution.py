"""The module for Genomic Substitution Validation."""
from typing import Optional, List
from .single_nucleotide_variation_base import SingleNucleotideVariationBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import GenomicSubstitutionToken
import logging
from variation.schemas.token_response_schema import Token

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicSubstitution(SingleNucleotideVariationBase):
    """The Genomic Substitution Validator class."""

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

    def get_hgvs_expr(self, classification, t, s, is_hgvs) -> str:
        """Return HGVS expression.

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :param Token s: The classification token
        :param bool is_hgvs: Whether or not classification is HGVS token
        :return: hgvs expression
        """
        if not is_hgvs:
            hgvs_expr = f"{t}:{s.reference_sequence.lower()}.{s.position}" \
                        f"{s.ref_nucleotide}>{s.new_nucleotide}"
        else:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type == 'HGVS'][0]
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
            'longest_compatible_remaining': dict(),
            'grch38': dict()
        }

        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                ref_nuc = \
                    self.seqrepo_access.sequence_at_position(t, s.position)

                allele, t, hgvs_expr, is_ensembl = \
                    self.get_allele_with_context(classification, t, s, errors)

                if allele:
                    self.check_ref_nucleotide(ref_nuc, s, t, errors)
                else:
                    errors.append("Unable to get allele")

                if not errors:
                    mane = self.mane_transcript.get_mane_transcript(
                        t, s.position, s.position, s.reference_sequence,
                        normalize_endpoint=normalize_endpoint
                    )

                    if mane:
                        if not gene_tokens and mane['gene']:
                            gene_tokens.append(
                                self._gene_matcher.match(mane['gene'])
                            )

                        if mane['status'] != 'grch38':
                            s.molecule_context = 'transcript'
                            s.reference_sequence = 'c'

                            # TODO: Only if on different strands
                            ref_rev = s.ref_nucleotide[::-1]
                            alt_rev = s.new_nucleotide[::-1]

                            complements = {
                                'A': 'T',
                                'T': 'A',
                                'C': 'G',
                                'G': 'C'
                            }

                            s.ref_nucleotide = ''
                            s.new_nucleotide = ''
                            for nt in ref_rev:
                                s.ref_nucleotide += complements[nt]
                            for nt in alt_rev:
                                s.new_nucleotide += complements[nt]

                        mane_hgvs_expr = \
                            f"{mane['refseq']}:{s.reference_sequence}." \
                            f"{mane['pos'][0]}{s.ref_nucleotide}>" \
                            f"{s.new_nucleotide}"

                        self.add_mane_data(
                            mane_hgvs_expr, mane, mane_data, s
                        )

                self.add_validation_result(
                    allele, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

                if is_hgvs:
                    break

        self.add_mane_to_validation_results(
            mane_data, valid_alleles, results, classification, gene_tokens
        )

    def get_gene_tokens(self, classification):
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
        return 'genomic substitution'

    def is_token_instance(self, t):
        """Check that token is genomic substitution."""
        return t.token_type == 'GenomicSubstitution'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is genomic
        substitution.
        """
        return classification_type == ClassificationType.GENOMIC_SUBSTITUTION

    def human_description(self, transcript,
                          token: GenomicSubstitutionToken) -> str:
        """Return a human description of the identified variation."""
        return f'A genomic substitution from {token.ref_nucleotide}' \
               f' to {token.new_nucleotide} at position ' \
               f'{token.position} on transcript {transcript}'
