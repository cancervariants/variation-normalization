"""The module for Polypeptide Sequence Variant Validation."""
from typing import List
from abc import abstractmethod
from .validator import Validator
from variant.schemas.classification_response_schema import Classification
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.validation_response_schema import ValidationResult, \
    LookupType
from variant.schemas.token_response_schema import Token
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import AminoAcidCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
import logging

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class PolypeptideSequenceVariantBase(Validator):
    """The Polypeptide Sequence Variant Validator Base class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 amino_acid_cache: AminoAcidCache) \
            -> None:
        """Initialize the validator.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        """
        super().__init__(seq_repo_access, transcript_mappings, gene_symbol)
        self._amino_acid_cache = amino_acid_cache

    def validate(self, classification: Classification) \
            -> List[ValidationResult]:
        """Validate a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :return: A list of validation results
        """
        results = list()
        errors = list()

        classification_tokens = [t for t in classification.all_tokens if self.is_token_instance(t)]  # noqa: E501
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

        transcripts = self.transcript_mappings.protein_transcripts(
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

    def get_hgvs_expr(self, classification) -> str:
        """Return HGVS expression for a classification.

        :param Classification classification: A classification for a list of
            tokens
        :return: The classification's HGVS expression as a string
        """
        # Replace `=` in silent mutation with 3 letter amino acid code
        hgvs_token = \
            [t for t in classification.all_tokens if
             isinstance(t, Token) and t.token_type == 'HGVS'][0]
        hgvs_expr = hgvs_token.input_string

        if '=' in hgvs_expr:
            hgvs_parsed = \
                self.hgvs_parser.parse_hgvs_variant(hgvs_expr)
            alt_amino = hgvs_parsed.posedit.pos.end.aa
            three_letter = \
                self._amino_acid_cache.amino_acid_code_conversion[alt_amino]
            hgvs_expr = hgvs_expr.replace('=', three_letter)
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

                if 'HGVS' in classification.matching_tokens:
                    hgvs_expr = self.get_hgvs_expr(classification)
                    allele = self.get_allele_from_hgvs(hgvs_expr, errors)
                    t = hgvs_expr.split(':')[0]
                else:
                    hgvs_expr = f"{t}:p.{s.ref_protein}{s.position}" \
                                f"{s.alt_protein}"
                    allele = self.get_allele_from_hgvs(hgvs_expr, errors)

                mane_transcripts_dict[hgvs_expr] = {
                    'classification_token': s,
                    'transcript_token': t
                }

                if allele and len(allele['state']['sequence']) == 3:
                    allele['state']['sequence'] = \
                        self._amino_acid_cache.convert_three_to_one(
                            allele['state']['sequence'])

                if not errors:
                    ref_protein = \
                        self.seqrepo_access.sequence_at_position(t, s.position)

                    if ref_protein and len(ref_protein) == 1 \
                            and len(s.ref_protein) == 3:
                        ref_protein = self._amino_acid_cache.amino_acid_code_conversion[ref_protein]  # noqa: E501
                    if ref_protein != s.ref_protein:
                        errors.append(f'Needed to find {s.ref_protein} at'
                                      f' position {s.position} on {t}'
                                      f' but found {ref_protein}')

                if not errors and allele not in valid_alleles:
                    results.append(self.get_validation_result(
                        classification, True, 1, allele,
                        self.human_description(t, s),
                        self.concise_description(t, s), [], gene_tokens))

                    valid_alleles.append(allele)
                else:
                    results.append(self.get_validation_result(
                        classification, False, 1, allele,
                        self.human_description(t, s),
                        self.concise_description(t, s), errors, gene_tokens))

        # Now add Mane transcripts to results
        self.add_mane_transcript(classification, results, gene_tokens,
                                 mane_transcripts_dict)

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_protein_gene_symbol_tokens(classification)

    def concise_description(self, transcript, token) -> str:
        """Return a description of the identified variant."""
        return f'{transcript} {token.ref_protein}' \
               f'{token.position}{token.alt_protein}'

    @abstractmethod
    def human_description(self, transcript, token) -> str:
        """Return a human description of the identified variant."""
        raise NotImplementedError
