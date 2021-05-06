"""The module for Amino Acid DelIns Validation."""
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import AminoAcidDelInsToken
from variant.schemas.validation_response_schema import LookupType
from typing import List
from variant.validators.validator import Validator
from variant.schemas.classification_response_schema import Classification
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.validation_response_schema import ValidationResult
from variant.schemas.token_response_schema import Token
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import AminoAcidCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
import logging


logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class AminoAcidDelIns(Validator):
    """The Amino Acid DelIns Validator class."""

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
                if 'HGVS' in classification.matching_tokens or \
                        'ReferenceSequence' in classification.matching_tokens:
                    hgvs_expr = self.get_hgvs_expr(classification, t, s, True)
                    allele = self.get_allele_from_hgvs(hgvs_expr, errors)
                    if allele:
                        t = hgvs_expr.split(':')[0]
                        # MANE Select Transcript for HGVS expressions
                        mane_transcripts_dict[hgvs_expr] = {
                            'classification_token': s,
                            'transcript_token': t
                        }
                else:
                    refseq = ([a for a in self.seqrepo_access.aliases(t)
                               if a.startswith('refseq:NP_')] or [None])[0]

                    if not refseq:
                        allele = None
                    else:
                        dels = f"{s.start_aa_del}{s.start_pos_del}"
                        if s.start_pos_del is not None and \
                                s.end_pos_del is not None:
                            dels += f"_{s.end_aa_del}{s.end_pos_del}"
                        hgvs_expr = f"{refseq.split('refseq:')[-1]}:p." \
                                    f"{dels}delins{s.inserted_sequence}"

                        # MANE Select Transcript for Gene Name + Variation
                        # (ex: BRAF V600E)
                        if hgvs_expr not in mane_transcripts_dict.keys():
                            mane_transcripts_dict[hgvs_expr] = {
                                'classification_token': s,
                                'transcript_token': t
                            }

                        allele = self.get_allele_from_transcript(
                            classification, t, s, errors
                        )

                if not allele:
                    errors.append("Unable to find allele.")
                else:
                    if len(allele['state']['sequence']) == 3:
                        allele['state']['sequence'] = \
                            self._amino_acid_cache.convert_three_to_one(
                                allele['state']['sequence'])

                # Check ref start/end protein matches expected
                if not errors:
                    ref_start_aa_del = \
                        self.seqrepo_access.sequence_at_position(t, s.start_pos_del)  # noqa: E501
                    if ref_start_aa_del and len(ref_start_aa_del) == 1 \
                            and len(s.start_aa_del):
                        ref_start_aa_del = \
                            self._amino_acid_cache.amino_acid_code_conversion[ref_start_aa_del]  # noqa: E501

                    if ref_start_aa_del != s.start_aa_del:
                        errors.append(f"Needed to find {s.start_aa_del} at "
                                      f"position {s.start_pos_del} on {t} "
                                      f"but found {ref_start_aa_del}")

                    if not errors and s.end_aa_del:
                        ref_end_aa_del = \
                            self.seqrepo_access.sequence_at_position(
                                t, s.end_pos_del
                            )

                        if ref_end_aa_del and len(ref_end_aa_del) == 1 \
                                and len(s.end_aa_del):
                            ref_end_aa_del = \
                                self._amino_acid_cache.amino_acid_code_conversion[ref_end_aa_del]  # noqa: E501

                        if ref_end_aa_del != s.end_aa_del:
                            errors.append(f"Needed to find {s.end_aa_del} "
                                          f"at position {s.end_pos_del} on "
                                          f"{t} but found {ref_end_aa_del}")

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

        # Now add MANE transcripts to result
        self.add_mane_transcript(classification, results, gene_tokens,
                                 mane_transcripts_dict)

    def get_allele_from_transcript(self, classification, t, s, errors):
        """Return allele from a given transcript.

        :param Classification classification: The classification for input str
        :param Classification s: Classification token
        :param str t: Transcript
        :param list errors: List of errors
        :return: Allele as a dictionary
        """
        hgvs_expr = self.get_hgvs_expr(classification, t, s, False)
        return self.get_allele_from_hgvs(hgvs_expr, errors)

    def get_hgvs_expr(self, classification, t, s, is_hgvs) -> str:
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
            dels = f"{s.start_aa_del}{s.start_pos_del}"
            if s.start_pos_del is not None and s.end_pos_del is not None:
                dels += f"_{s.end_aa_del}{s.end_pos_del}"
            hgvs_expr = f"{prefix}{dels}delins{s.inserted_sequence}"
        else:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type == 'HGVS']
            if not hgvs_token:
                hgvs_token = \
                    [t for t in classification.all_tokens if
                     isinstance(t, Token) and t.token_type == 'ReferenceSequence']  # noqa: E501
            hgvs_token = hgvs_token[0]
            hgvs_expr = hgvs_token.input_string
        return hgvs_expr

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_protein_gene_symbol_tokens(classification)

    def variant_name(self):
        """Return the variant name."""
        return 'amino acid delins'

    def is_token_instance(self, t):
        """Check that token is Amino Acid DelIns."""
        return t.token_type == 'AminoAcidDelIns'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Amino Acid DelIns.
        """
        return classification_type == ClassificationType.AMINO_ACID_DELINS

    def concise_description(self, transcript, token) -> str:
        """Return a description of the identified variant."""
        dels = f"{token.start_aa_del}{token.start_pos_del}"
        if token.start_pos_del is not None and token.end_pos_del is not None:
            dels += f"_{token.end_aa_del}{token.end_pos_del}"

        return f'{transcript}:p.{dels}delins{token.inserted_sequence}'

    def human_description(self, transcript,
                          token: AminoAcidDelInsToken) -> str:
        """Return a human description of the identified variant."""
        position = f"{token.start_aa_del}{token.start_pos_del}"
        if token.start_pos_del is not None and token.end_pos_del is not None:
            position += f"to {token.end_aa_del}{token.end_pos_del}"

        return f"An Amino Acid DelIns deletion of {position} replaced by " \
               f"{token.input_string} on transcript {transcript}"
