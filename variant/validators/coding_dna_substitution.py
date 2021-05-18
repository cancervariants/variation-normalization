"""The module for Coding DNA Substitution Validation."""
from .single_nucleotide_variant_base import SingleNucleotideVariantBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import CodingDNASubstitutionToken
from typing import List, Optional
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.token_response_schema import Token
import logging

# TODO:
#  LRG_ (LRG_199t1:c)


logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class CodingDNASubstitution(SingleNucleotideVariantBase):
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
                          isinstance(t, Token) and t.token_type == 'HGVS'][0]
            hgvs_expr = hgvs_token.input_string
        return hgvs_expr

    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results, gene_tokens) \
            -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of Tokens
        :param list transcripts: A list of transcript strings
        :param Classification classification: A classification for a list of
            tokens
        :param list results: A list to store validation result objects
        :param list gene_tokens: List of GeneMatchTokens
        """
        valid_alleles = list()
        mane_transcripts_dict = dict()
        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                allele, t, hgvs_expr, is_ensembl = \
                    self.get_allele_with_context(classification, t, s, errors)

                if not allele:
                    errors.append("Unable to find allele.")
                else:
                    mane_transcripts_dict[hgvs_expr] = {
                        'classification_token': s,
                        'transcript_token': t,
                        'nucleotide': is_ensembl
                    }

                    ref_nuc = \
                        self.seqrepo_access.sequence_at_position(t, s.position)
                    self.check_ref_nucleotide(ref_nuc, s, t, errors)

                self.add_validation_result(
                    allele, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

        # Now add Mane transcripts to results
        self.add_mane_transcript(classification, results, gene_tokens,
                                 mane_transcripts_dict)

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_coding_dna_gene_symbol_tokens(classification)

    def variant_name(self):
        """Return the variant name."""
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
        """Return a human description of the identified variant."""
        return f'A coding DNA substitution from {token.ref_nucleotide}' \
               f' to {token.new_nucleotide} at position ' \
               f'{token.position} on transcript {transcript}'
