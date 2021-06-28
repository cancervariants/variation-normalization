"""The module for Genomic DelIns Validation."""
from variation.validators.delins_base import DelInsBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import GenomicDelInsToken
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.token_response_schema import Token
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicDelIns(DelInsBase):
    """The Genomic DelIns Validator class."""

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
        """Return HGVS expression

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :param Token s: The classification token
        :param bool is_hgvs: Whether or not classification is HGVS token
        :return: hgvs expression
        """
        if not is_hgvs:
            prefix = f"{t}:{s.reference_sequence.lower()}."
            if s.start_pos_del is not None and s.end_pos_del is not None:
                pos_del = f"{s.start_pos_del}_{s.end_pos_del}"
            else:
                pos_del = s.start_pos_del

            if s.inserted_sequence1 is not None and \
                    s.inserted_sequence2 is not None:
                inserted_seq = f"{s.inserted_sequence1}_{s.inserted_sequence2}"
            else:
                inserted_seq = s.inserted_sequence1

            hgvs_expr = f"{prefix}{pos_del}delins{inserted_seq}"
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
            'longest_compatible_remaining': dict()
        }

        for s in classification_tokens:
            for t in transcripts:
                errors = list()

                allele, t, hgvs_expr, is_ensembl = \
                    self.get_allele_with_context(classification, t, s, errors)

                mane = self.mane_transcript.get_mane_transcript(
                    t, s.start_pos_del, s.end_pos_del, s.reference_sequence,
                    normalize_endpoint=normalize_endpoint
                )
                # TODO: Fix MANE when GRCh38 rather than mane
                if mane:
                    if not gene_tokens:
                        gene_tokens.append(
                            self._gene_matcher.match(mane['gene'])
                        )

                    prefix = f"{mane['refseq']}:c."
                    if s.start_pos_del is not None and s.end_pos_del is not None:  # noqa: E501
                        pos_del = f"{mane['pos'][0]}_{mane['pos'][1]}"
                    else:
                        pos_del = f"{mane['pos'][0]}"
                    if s.inserted_sequence1 is not None \
                            and s.inserted_sequence2 is not None:
                        inserted_seq = \
                            f"{s.inserted_sequence1}_{s.inserted_sequence2}"
                    else:
                        inserted_seq = f"{s.inserted_sequence1}"
                    mane_hgvs_expr = f"{prefix}{pos_del}delins{inserted_seq}"
                    self.add_mane_data(mane_hgvs_expr, mane, mane_data, s)

                if allele:
                    self.check_pos_index(t, s, errors)

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
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
        return 'genomic delins'

    def is_token_instance(self, t):
        """Check that token is Genomic DelIns."""
        return t.token_type == 'GenomicDelIns'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Genomic DelIns.
        """
        return classification_type == ClassificationType.GENOMIC_DELINS

    def human_description(self, transcript,
                          token: GenomicDelInsToken) -> str:
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

        return f"A Genomic DelIns deletion of {position} replaced by " \
               f"{sequence} on transcript {transcript}"
