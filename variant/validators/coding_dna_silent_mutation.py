"""The module for Coding DNA Substitution Validation."""
from .single_nucleotide_variant_base import SingleNucleotideVariantBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import CodingDNASilentMutationToken
from variant.schemas.validation_response_schema import LookupType
from typing import List
from variant.schemas.classification_response_schema import Classification
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.validation_response_schema import ValidationResult
from variant.schemas.token_response_schema import Token
import logging


logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class CodingDNASilentMutation(SingleNucleotideVariantBase):
    """The Coding DNA Silent Mutation Validator class."""

    def validate(self, classification: Classification) \
            -> List[ValidationResult]:
        """Validate a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :return: A list of validation results
        """
        results = list()
        errors = list()

        classification_tokens = self.get_classification_tokens(classification)
        gene_tokens = self.get_gene_tokens(classification)

        if len(classification.non_matching_tokens) > 0:
            errors.append(f"Non matching tokens found for "
                          f"{self.variant_name()}.")

        if len(gene_tokens) == 0:
            errors.append(f'No gene tokens for a {self.variant_name()}.')

        if len(gene_tokens) > 1:
            errors.append('More than one gene symbol found for a single'
                          f' {self.variant_name()}')

        if len(errors) > 0:
            return [self.get_validation_result(
                classification, False, 0, None,
                '', '', errors, gene_tokens)]

        transcripts = self.transcript_mappings.coding_dna_transcripts(
            gene_tokens[0].token, LookupType.GENE_SYMBOL)

        if not transcripts:
            errors.append(f'No transcripts found for gene symbol '
                          f'{gene_tokens[0].token}')
            return [self.get_validation_result(
                classification, False, 0, None,
                '', '', errors, gene_tokens)]

        self.get_valid_invalid_results(classification_tokens, transcripts,
                                       classification, results, gene_tokens)
        return results

    def get_hgvs_expr(self, classification, t, s, is_hgvs):
        """Return HGVS expression and whether or not it's an Ensembl transcript

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
         :param bool is_hgvs: Whether or not classification is HGVS token
        :return: A tuple containing the hgvs expression and whether or not
            it's an Ensembl Transcript
        """
        # Get transcript
        if is_hgvs:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type == 'HGVS'][0]
            input_string = hgvs_token.input_string
            if not input_string.startswith('ENST'):
                t = input_string.split(':')[0]

        # Make hgvs_expr
        ref_nuc = self.seqrepo_access.sequence_at_position(t, s.position)
        if not ref_nuc:
            hgvs_expr = None
        else:
            nucleotides = {'C', 'T', 'A', 'G'}
            new_nuc = (nucleotides - {ref_nuc}).pop()
            hgvs_expr = f"{t}:{s.reference_sequence.lower()}." \
                        f"{s.position}{ref_nuc}>{new_nuc}"

        return hgvs_expr, None

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
        for s in classification_tokens:
            for t in transcripts:
                errors = list()

                if 'HGVS' in classification.matching_tokens:
                    # TODO: How to convert ENST_ to NM_ versioned
                    hgvs_expr, _ = self.get_hgvs_expr(classification,
                                                      t, s, True)
                    allele = self.get_allele_from_hgvs(hgvs_expr, errors)
                    if hgvs_expr:
                        t = hgvs_expr.split(':')[0]
                else:
                    allele = self.get_allele_from_transcript(classification,
                                                             t, s, errors)

                if allele:
                    # Fix ref nuc in allele
                    ref_nuc = \
                        self.seqrepo_access.sequence_at_position(t, s.position)
                    allele['state']['sequence'] = ref_nuc

                    len_of_seq = self.seqrepo_access.len_of_sequence(t)
                    if len_of_seq < s.position - 1:
                        errors.append('Sequence index error')

                self.add_validation_result(
                    allele, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_coding_dna_gene_symbol_tokens(classification)

    def variant_name(self):
        """Return the variant name."""
        return 'coding dna silent mutation'

    def is_token_instance(self, t):
        """Check that token is Coding DNA Silent Mutation."""
        return t.token_type == 'CodingDNASilentMutation'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is coding dna silent
        mutation.
        """
        return classification_type == ClassificationType.CODING_DNA_SILENT_MUTATION  # noqa: E501

    def human_description(self, transcript,
                          token: CodingDNASilentMutationToken) -> str:
        """Return a human description of the identified variant."""
        return f'A coding DNA silent mutation from {token.position} ' \
               'was a {token.ref_nucleotide} (the nucleotide was not changed)'
