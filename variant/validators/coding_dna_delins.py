"""The module for Coding DNA DelIns Validation."""
from variant.validators.delins_base import DelInsBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import CodingDNADelInsToken
from variant.schemas.validation_response_schema import LookupType
from typing import List, Tuple
from variant.schemas.classification_response_schema import Classification
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.validation_response_schema import ValidationResult
from variant.schemas.token_response_schema import Token
import logging


logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class CodingDNADelIns(DelInsBase):
    """The Coding DNA DelIns Validator class."""

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

    def get_hgvs_expr(self, classification, t, s, is_hgvs) -> Tuple[str, bool]:
        """Return HGVS expression and whether or not it's an Ensembl transcript

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :param bool is_hgvs: Whether or not classification is HGVS token
        :return: A tuple containing the hgvs expression and whether or not
            it's an Ensembl Transcript
        """
        if not is_hgvs:
            prefix = f"{t}:{s.reference_sequence.lower()}."
            if s.start_pos_del is not None and s.end_pos_del is not None:
                pos_del = f"{s.start_pos_del}_{s.end_pos_del}"
            else:
                pos_del = s.end_pos_del

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

        gene_token = [t for t in classification.all_tokens
                      if t.token_type == 'GeneSymbol']
        if gene_token:
            is_ensembl_transcript = True
        else:
            is_ensembl_transcript = False
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

                if 'HGVS' in classification.matching_tokens:
                    hgvs_expr, is_ensembl_transcript = \
                        self.get_hgvs_expr(classification, t, s, True)
                    allele = self.get_allele_from_hgvs(hgvs_expr, errors)
                    if allele:
                        t = hgvs_expr.split(':')[0]
                        # MANE Select Transcript for HGVS expressions
                        mane_transcripts_dict[hgvs_expr] = {
                            'classification_token': s,
                            'transcript_token': t,
                            'is_ensembl_transcript': is_ensembl_transcript
                        }
                    else:
                        # Typically cant find an allele when querying ENST
                        # So try again using transcript mapping
                        errors = list()
                        allele = self.get_allele_from_transcript(
                            classification, t, s, errors
                        )

                        if allele:
                            # MANE Select Transcript for Gene Name + Variation
                            # (ex: BRAF V600E)
                            self._add_hgvs_to_mane_transcripts_dict(
                                classification, mane_transcripts_dict, s, t,
                                gene_tokens
                            )
                else:
                    allele = self.get_allele_from_transcript(classification,
                                                             t, s, errors)
                    if allele:
                        # MANE Select Transcript for Gene Name + Variation
                        # (ex: BRAF V600E)
                        self._add_hgvs_to_mane_transcripts_dict(
                            classification, mane_transcripts_dict, s, t,
                            gene_tokens
                        )

                if allele:
                    len_of_seq = self.seqrepo_access.len_of_sequence(t)
                    is_len_lt_end = len_of_seq < int(s.end_pos_del) - 1
                    is_len_lt_start = \
                        s.start_pos_del and len_of_seq < int(s.start_pos_del) - 1  # noqa: E501

                    if is_len_lt_end or is_len_lt_start:
                        errors.append('Sequence index error')

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
        """Return a human description of the identified variant."""
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
