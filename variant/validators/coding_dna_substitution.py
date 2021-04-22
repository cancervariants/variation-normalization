"""The module for Coding DNA Substitution Validation."""
from .single_nucleotide_variant_base import SingleNucleotideVariantBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import CodingDNASubstitutionToken
from variant.schemas.validation_response_schema import LookupType
from typing import List
from variant.schemas.classification_response_schema import Classification
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.validation_response_schema import ValidationResult
from variant.schemas.token_response_schema import Token
import logging

# TODO:
#  LRG_ (LRG_199t1:c)


logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class CodingDNASubstitution(SingleNucleotideVariantBase):
    """The Coding DNA Substitution Validator class."""

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

    def get_hgvs_expr(self, classification, t) -> tuple:
        """Return HGVS expression and whether or not it's an Ensembl transcript

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :return: A tuple containing the hgvs expression and whether or not
            it's an Ensembl Transcript
        """
        hgvs_token = \
            [t for t in classification.all_tokens if
             isinstance(t, Token) and t.token_type == 'HGVS'][0]
        hgvs_expr = hgvs_token.input_string
        if hgvs_expr.startswith('ENST'):
            is_ensembl_transcript = True
            if not t.startswith('ENST'):
                hgvs_expr = f"{t}:{hgvs_expr.split(':')[1]}"
        else:
            is_ensembl_transcript = False

        gene_token = [t for t in classification.all_tokens
                      if t.token_type == 'GeneSymbol']
        if gene_token:
            is_ensembl_transcript = True

        # Replace `=` in silent mutation
        if '=' in hgvs_expr:
            # TODO
            pass
        return hgvs_expr, is_ensembl_transcript

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
                ref_nuc = \
                    self.seqrepo_access.sequence_at_position(t, s.position)

                if 'HGVS' in classification.matching_tokens and \
                        not t.startswith('ENST'):
                    hgvs_expr, is_ensembl_transcript = \
                        self.get_hgvs_expr(classification, t)
                    allele = self.get_allele_from_hgvs(hgvs_expr, errors)
                    if allele:
                        # MANE Select Transcript for HGVS expressions
                        mane_transcripts_dict[hgvs_expr] = {
                            'classification_token': s,
                            'transcript_token': t,
                            'is_ensembl_transcript': is_ensembl_transcript
                        }
                else:
                    allele = self.get_allele_from_transcript(s, t, errors)
                    if allele:
                        # MANE Select Transcript for Gene Name + Variation
                        # (ex: BRAF V600E)
                        self._add_hgvs_to_mane_transcripts_dict(
                            classification, mane_transcripts_dict, s, t,
                            gene_tokens
                        )

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
        gene_tokens = self.get_gene_symbol_tokens(classification)
        if not gene_tokens:
            # Convert refseq to gene symbol
            refseq = \
                ([t.token for t in classification.all_tokens if
                 t.token_type in ['HGVS', 'ReferenceSequence']] or [None])[0]

            if not refseq:
                return []

            if ':' in refseq:
                refseq = refseq.split(':')[0]

            res = self.seqrepo_access.aliases(refseq)
            aliases = [a.split('refseq:')[1] for a
                       in res if a.startswith('refseq')]

            if not aliases:
                aliases = [refseq]

            gene_symbols = list()
            if aliases:
                for alias in aliases:
                    gene_symbol = \
                        self.transcript_mappings.get_gene_symbol_from_refseq_rna(alias)  # noqa: E501

                    if not gene_symbol:
                        gene_symbol = \
                            self.transcript_mappings.get_gene_symbol_from_ensembl_transcript(alias)  # noqa: E501

                    if gene_symbol:
                        if gene_symbol not in gene_symbols:
                            gene_symbols.append(gene_symbol)
                            gene_tokens.append(
                                self._gene_matcher.match(gene_symbol)
                            )
                    else:
                        logger.warning(f"No gene symbol found for rna "
                                       f"{alias} in transcript_mappings.tsv")
        return gene_tokens

    def variant_name(self):
        """Return the variant name."""
        return 'coding dna substitution'

    def is_token_instance(self, t):
        """Check that token is Coding DNA Substitution."""
        return t.token_type == 'CodingDNASubstitution'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is amino acid
        substitution.
        """
        return classification_type == ClassificationType.CODING_DNA_SUBSTITUTION  # noqa: E501

    def human_description(self, transcript,
                          psub_token: CodingDNASubstitutionToken) -> str:
        """Return a human description of the identified variant."""
        return f'An coding DNA substitution from {psub_token.ref_nucleotide}' \
               f' to {psub_token.new_nucleotide} at position ' \
               f'{psub_token.position} on transcript {transcript}'
