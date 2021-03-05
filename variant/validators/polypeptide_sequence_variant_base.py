"""The module for Polypeptide Sequence Variant Validation."""
from typing import List
from abc import abstractmethod
from requests.exceptions import HTTPError
from .validator import Validator
from variant.schemas.classification_response_schema import Classification
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.validation_response_schema import ValidationResult, \
    LookupType
from variant.schemas.token_response_schema import Token
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import AminoAcidCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
import hgvs.parser


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
        self.transcript_mappings = transcript_mappings
        self.seq_repo_access = seq_repo_access
        self._gene_matcher = gene_symbol
        self._amino_acid_cache = amino_acid_cache
        self.dp = SeqRepoDataProxy(seq_repo_access.seq_repo_client)
        self.tlr = Translator(data_proxy=self.dp)
        self.hgvs_parser = hgvs.parser.Parser()

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

    def get_vrs_allele(self, sequence_id, s) -> dict:
        """Return VRS Allele object.

        :param str sequence_id: Sequence containing the sequence to be located
        :param Token s: A token
        :return: A VRS Allele object as a dictionary
        """
        seq_location = models.SequenceLocation(
            sequence_id=sequence_id,
            interval=models.SimpleInterval(
                start=s.position - 1,
                end=s.position
            )
        )

        state = models.SequenceState(sequence=s.alt_protein.upper())
        state_dict = state.as_dict()
        if len(state_dict['sequence']) == 3:
            state.sequence = self._amino_acid_cache._convert_three_to_one(
                state_dict['sequence'])

        allele = models.Allele(location=seq_location,
                               state=state)
        allele['_id'] = ga4gh_identify(allele)
        return allele.as_dict()

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
                self._amino_acid_cache._amino_acid_code_conversion[alt_amino]
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
        for s in classification_tokens:
            for t in transcripts:
                valid = True
                errors = list()
                ref_protein = \
                    self.seq_repo_access.protein_at_position(t, s.position)

                if 'HGVS' in classification.matching_tokens:
                    hgvs_expr = self.get_hgvs_expr(classification)
                    try:
                        # TODO: We should get this to work w/o doing above
                        #       replacement
                        seq_location = self.tlr.translate_from(hgvs_expr,
                                                               'hgvs')
                    except HTTPError:
                        errors.append(f"{hgvs_expr} is an invalid "
                                      f"HGVS expression.")
                        valid = False
                        allele = None
                    else:
                        allele = seq_location.as_dict()
                else:
                    sequence_id = \
                        self.dp.translate_sequence_identifier(t, 'ga4gh')[0]
                    allele = self.get_vrs_allele(sequence_id, s)

                if len(allele['state']['sequence']) == 3:
                    allele['state']['sequence'] = \
                        self._amino_acid_cache._convert_three_to_one(
                            allele['state']['sequence'])

                if not errors:
                    if ref_protein and len(ref_protein) == 1 \
                            and len(s.ref_protein) == 3:
                        ref_protein = self._amino_acid_cache._amino_acid_code_conversion[ref_protein]  # noqa: E501
                    if ref_protein != s.ref_protein:
                        errors.append(f'Needed to find {s.ref_protein} at'
                                      f' position {s.position} on {t}'
                                      f' but found {ref_protein}')
                        valid = False

                if valid and allele not in valid_alleles:
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

    def get_validation_result(self, classification, is_valid, confidence_score,
                              allele, human_description, concise_description,
                              errors, gene_tokens) -> ValidationResult:
        """Return a validation result object.

        :param Classification classification: The classification for tokens
        :param boolean is_valid: Whether or not the classification is valid
        :param int confidence_score: The classification confidence score
        :param dict allele: A VRS Allele object
        :param str human_description: A human description describing the
            variant
        :param str concise_description: The identified variant
        :param list errors: A list of errors for the classification
        :param list gene_tokens: List of GeneMatchTokens
        :return: A validation result
        """
        return ValidationResult(
            classification=classification,
            is_valid=is_valid,
            confidence_score=confidence_score,
            allele=allele,
            human_description=human_description,
            concise_description=concise_description,
            errors=errors,
            gene_tokens=gene_tokens
        )

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        gene_tokens = [t for t in classification.all_tokens
                       if t.token_type == 'GeneSymbol']
        if not gene_tokens:
            # Convert refseq to gene symbol
            refseq = \
                [t.token for t in classification.all_tokens if
                 t.token_type in ['HGVS', 'ReferenceSequence']][0]

            if ':' in refseq:
                refseq = refseq.split(':')[0]

            res = self.seq_repo_access.aliases(refseq)
            aliases = [a.split('ensembl:')[1] for a
                       in res if a.startswith('ensembl')]

            gene_symbols = list()
            for alias in aliases:
                gene_symbol = \
                    self.transcript_mappings.refseq_gene_symbol(alias)

                if gene_symbol and gene_symbol not in gene_symbols:
                    gene_symbols.append(gene_symbol)
                    gene_tokens.append(self._gene_matcher.match(gene_symbol))
        return gene_tokens

    @abstractmethod
    def is_token_instance(self, t) -> bool:
        """Check to see if token is instance of a token type."""
        raise NotImplementedError

    @abstractmethod
    def variant_name(self):
        """Return the variant name"""
        raise NotImplementedError

    @abstractmethod
    def validates_classification_type(self, classification_type) \
            -> bool:
        """Return the classification type"""
        raise NotImplementedError

    def concise_description(self, transcript, token) -> str:
        """Return a description of the identified variant."""
        return f'{transcript} {token.ref_protein}' \
               f'{token.position}{token.alt_protein}'

    @abstractmethod
    def human_description(self, transcript, token) -> str:
        """Return a human description of the identified variant."""
        raise NotImplementedError
