"""Module for Validation."""
import re
from typing import List, Tuple, Optional
from abc import ABC, abstractmethod
from variation.schemas.classification_response_schema import Classification, \
    ClassificationType
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.validation_response_schema import ValidationResult, \
    LookupType
from variation.tokenizers import GeneSymbol
from variation.data_sources import SeqRepoAccess, TranscriptMappings
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
import hgvs.parser
from requests.exceptions import HTTPError
import logging
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from variation.validators.genomic_base import GenomicBase
from variation.data_sources import UTA

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class Validator(ABC):
    """The validator class."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTA) -> None:
        """Initialize the DelIns validator.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: Gene symbol tokenizer
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param UTA uta: Access to UTA queries
        """
        self.transcript_mappings = transcript_mappings
        self.seqrepo_access = seqrepo_access
        self._gene_matcher = gene_symbol
        self.dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        self.tlr = Translator(data_proxy=self.dp)
        self.hgvs_parser = hgvs.parser.Parser()
        self.uta = uta
        self.genomic_base = GenomicBase(self.dp, self.uta)
        self.mane_transcript = mane_transcript

    @abstractmethod
    def is_token_instance(self, t) -> bool:
        """Check to see if token is instance of a token type.

        :param Token t: Classification token to find type of
        :return: `True` if token is instance of class token. `False` otherwise.
        """
        raise NotImplementedError

    @abstractmethod
    def variation_name(self) -> str:
        """Return the variation name.

        :return: variation class name
        """
        raise NotImplementedError

    @abstractmethod
    def human_description(self, transcript, token) -> str:
        """Return a human description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: Human description of the variation change
        """
        raise NotImplementedError

    @abstractmethod
    def concise_description(self, transcript, token) -> str:
        """Return a HGVS description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression
        """
        raise NotImplementedError

    @abstractmethod
    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return a list of gene tokens for a classification.

        :param Classification classification: Classification for a list of
            tokens
        :return: A list of gene tokens for the classification
        """
        raise NotImplementedError

    @abstractmethod
    def get_hgvs_expr(self, classification, t, s, is_hgvs) -> str:
        """Return HGVS expression for a classification token.

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript accession
        :param Token s: The classification token
        :param bool is_hgvs: Whether or not classification is a HGVS expression
            token
        :return: Tuple[hgvs expression, `True` if MANE transcript should use
            Ensembl accession. `False` if MANE transcript should use RefSeq
            accession
        """
        raise NotImplementedError

    @abstractmethod
    def get_transcripts(self, gene_tokens, classification, errors)\
            -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param list gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param list errors: List of errors
        :return: List of transcript accessions
        """
        raise NotImplementedError

    @abstractmethod
    def validates_classification_type(self,
                                      classification_type: ClassificationType)\
            -> bool:
        """Check that classification type can be validated by validator.

        :param str classification_type: The type of variation
        :return: `True` if classification_type matches validator's
            classification type. `False` otherwise.
        """
        raise NotImplementedError

    @abstractmethod
    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results, gene_tokens,
                                  normalize_endpoint) -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of classification Tokens
        :param list transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param list results: Stores validation result objects
        :param list gene_tokens: List of GeneMatchTokens for a classification
        """
        raise NotImplementedError

    def validate(self, classification: Classification, normalize_endpoint) \
            -> List[ValidationResult]:
        """Return validation result for a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :return: Validation Result's containing valid and invalid results
        """
        results = list()
        errors = list()

        classification_tokens = self.get_classification_tokens(classification)
        if len(classification.non_matching_tokens) > 0:
            errors.append(f"Non matching tokens found for "
                          f"{self.variation_name()}.")

        gene_tokens = self.get_gene_tokens(classification)
        if len(gene_tokens) > 1:
            errors.append('More than one gene symbol found for a single'
                          f' {self.variation_name()}')

        try:
            # NC_ queries do not have gene tokens
            transcripts = \
                self.get_transcripts(gene_tokens, classification, errors)
        except IndexError:
            transcripts = list()

        if len(errors) > 0:
            return [self.get_validation_result(
                classification, False, 0, None,
                '', '', errors, gene_tokens)]

        self.get_valid_invalid_results(
            classification_tokens, transcripts, classification,
            results, gene_tokens, normalize_endpoint
        )
        return results

    def get_validation_result(self, classification, is_valid, confidence_score,
                              allele, human_description, concise_description,
                              errors, gene_tokens,
                              mane_transcript=None) -> ValidationResult:
        """Return a validation result object.

        :param Classification classification: The classification for tokens
        :param boolean is_valid: Whether or not the classification is valid
        :param int confidence_score: The classification confidence score
        :param dict allele: A VRS Allele object
        :param str human_description: A human description describing the
            variation
        :param str concise_description: HGVS expression for variation
        :param list errors: A list of errors for the classification
        :param list gene_tokens: List of GeneMatchTokens
        :param str mane_transcript: HGVS expression for MANE transcript found
            from ClinGen Allele Registry API
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
            gene_tokens=gene_tokens,
            mane_transcript=mane_transcript
        )

    def get_protein_transcripts(self, gene_tokens, errors)\
            -> Optional[List[str]]:
        """Get transcripts for variations with protein reference sequence.

        :param list gene_tokens: List of gene tokens for a classification
        :param list errors: List of errors
        :return: List of possible transcript accessions for the variation
        """
        transcripts = self.transcript_mappings.protein_transcripts(
            gene_tokens[0].token, LookupType.GENE_SYMBOL)
        if not transcripts:
            errors.append(f'No transcripts found for gene symbol '
                          f'{gene_tokens[0].token}')
        return transcripts

    def get_coding_dna_transcripts(self, gene_tokens, errors)\
            -> Optional[List[str]]:
        """Get transcripts for variations with coding DNA reference sequence.

        :param list gene_tokens: List of gene tokens for a classification
        :param list errors: List of errors
        :return: List of possible transcript accessions for the variation
        """
        transcripts = self.transcript_mappings.coding_dna_transcripts(
            gene_tokens[0].token, LookupType.GENE_SYMBOL)
        if not transcripts:
            errors.append(f'No transcripts found for gene symbol '
                          f'{gene_tokens[0].token}')
        return transcripts

    def get_genomic_transcripts(self, classification, errors)\
            -> Optional[List[str]]:
        """Get NC accessions for variations with genomic reference sequence.

        :param Classification classification: Classification for a list of
            tokens
        :param list errors: List of errors
        :return: List of possible NC accessions for the variation
        """
        nc_accessions = self.genomic_base.get_nc_accessions(classification)
        if not nc_accessions:
            errors.append('Could not find NC_ accession for '
                          f'{self.variation_name()}')
        return nc_accessions

    def get_classification_tokens(self, classification)\
            -> Optional[List[Classification]]:
        """Get classification tokens for a given instance.

        :param Classification classification: A classification for a list of
            tokens
        :return: A list of classification tokens
        """
        return [t for t in classification.all_tokens
                if self.is_token_instance(t)]

    def get_gene_symbol_tokens(self, classification)\
            -> Optional[List[GeneMatchToken]]:
        """Return tokens with GeneSymbol token type from a classification.

        :param Classification classification: Classification of input string
        :return: List of Gene Match Tokens
        """
        return [t for t in classification.all_tokens
                if t.token_type == 'GeneSymbol']

    def _add_gene_symbol_to_tokens(self, gene_symbol, gene_symbols,
                                   gene_tokens) -> None:
        """Add a gene symbol to list of gene match tokens.

        :param str gene_symbol: Gene symbol
        :param list gene_symbols: List of gene symbols matched
        :param list gene_tokens: List of GeneMatchTokens
        """
        if gene_symbol and gene_symbol not in gene_symbols:
            gene_symbols.append(gene_symbol)
            gene_tokens.append(self._gene_matcher.match(
                gene_symbol))

    def _get_gene_tokens(self, classification, mappings)\
            -> Optional[List[GeneMatchToken]]:
        """Get gene symbol tokens for protein or transcript reference
        sequences.

        :param Classification classification: Classification for a list of
            tokens
        :param list mappings: List of transcript mapping methods for
            corresponding reference sequence
        :return: A list of gene match tokens
        """
        gene_tokens = self.get_gene_symbol_tokens(classification)
        if not gene_tokens:
            refseq = \
                ([t.token for t in classification.all_tokens if
                  t.token_type in ['HGVS', 'ReferenceSequence',
                                   'LocusReferenceGenomic']] or [None])[0]

            if not refseq:
                return []

            if ':' in refseq:
                refseq = refseq.split(':')[0]

            gene_symbols = list()
            for mapping in mappings:
                gene_symbol = mapping(refseq)
                self._add_gene_symbol_to_tokens(
                    gene_symbol, gene_symbols, gene_tokens
                )
                if gene_tokens:
                    break
        return gene_tokens

    def get_protein_gene_symbol_tokens(self, classification)\
            -> Optional[List[GeneMatchToken]]:
        """Return gene tokens for a classification with protein reference
        sequence.

        :param Classification classification: The classification for a list of
            tokens
        :return: A list of Gene Match Tokens in the classification
        """
        mappings = [
            self.transcript_mappings.get_gene_symbol_from_ensembl_protein,
            self.transcript_mappings.get_gene_symbol_from_refeq_protein,
            self.transcript_mappings.get_gene_symbol_from_lrg
        ]
        return self._get_gene_tokens(classification, mappings)

    def get_coding_dna_gene_symbol_tokens(self, classification)\
            -> Optional[List[GeneMatchToken]]:
        """Return gene symbol tokens for classifications with coding dna
        reference sequence.

        :param Classification classification: Classification of input string
        :return: A list of gene match tokens
        """
        mappings = [
            self.transcript_mappings.get_gene_symbol_from_refseq_rna,
            self.transcript_mappings.get_gene_symbol_from_ensembl_transcript,  # noqa: E501
            self.transcript_mappings.get_gene_symbol_from_lrg
        ]
        return self._get_gene_tokens(classification, mappings)

    def get_allele_with_context(self, classification, t, s, errors)\
            -> Tuple[dict, str, str, bool]:
        """Get VRS allele object, transcript accession used, HGVS expression
        for transcript, and whether or not to use Ensembl MANE transcript.

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript accession retrieved from transcript mapping
        :param Token s: The classification token
        :param list errors: List of errors
        :return: Tuple[VRS allele, transcript accession, hgvs_expression,
            `True` if MANE should use Ensembl accession.
            `False` if MANE should use RefSeq accession]
        """
        if 'HGVS' in classification.matching_tokens or \
                'ReferenceSequence' in classification.matching_tokens:
            hgvs_expr = self.get_hgvs_expr(classification, t, s, True)
            t = hgvs_expr.split(':')[0]
        else:
            hgvs_expr = self.get_hgvs_expr(classification, t, s, False)
        allele = self.get_allele_from_hgvs(hgvs_expr, errors)

        gene_token = [t for t in classification.all_tokens
                      if t.token_type == 'GeneSymbol']
        if gene_token or hgvs_expr.startswith('EN'):
            is_ensembl = True
        else:
            is_ensembl = False
        return allele, t, hgvs_expr, is_ensembl

    def get_allele_from_hgvs(self, hgvs_expr, errors) -> Optional[dict]:
        """Return allele given a HGVS expression.

        :param str hgvs_expr: The HGVS expression
        :param list errors: List of errors
        :return: VRS Allele object represented as a dictionary
        """
        allele = None
        try:
            allele = self.tlr.translate_from(hgvs_expr, 'hgvs')
        except HTTPError:
            errors.append(f"{hgvs_expr} is an invalid HGVS expression.")
        except KeyError:
            errors.append("GA4GH Data Proxy unable to translate sequence "
                          f"identifier: {hgvs_expr}.")
        except ValueError:
            errors.append(f"Unable to parse {hgvs_expr} as hgvs variation")
        except:  # noqa
            errors.append("Unable to get VRS Allele.")
        else:
            allele = allele.as_dict()
            allele_seq_id = allele['location']['sequence_id']
            if allele_seq_id.startswith('ga4gh:GS.'):
                allele['location']['sequence_id'] = \
                    allele_seq_id.replace('ga4gh:GS.', 'ga4gh:SQ.')
        return allele

    def add_validation_result(self, allele, valid_alleles, results,
                              classification, s, t, gene_tokens, errors,
                              mane_transcript=None) -> bool:
        """Add validation result to list of results.

        :param dict allele: A VRS Allele object
        :param list valid_alleles: A list containing current valid alleles
        :param list results: A list of validation results
        :param Classification classification: The classification for tokens
        :param Token s: The classification token
        :param string t: Transcript
        :param list gene_tokens: List of GeneMatchTokens
        :param list errors: A list of errors for the classification
        :param str mane_transcript: The mane transcript found
        """
        if not errors:
            if mane_transcript or (allele and allele not in valid_alleles):
                results.append(
                    self.get_validation_result(
                        classification, True, 1, allele,
                        self.human_description(t, s),
                        self.concise_description(t, s), [],
                        gene_tokens, mane_transcript
                    )
                )
                valid_alleles.append(allele)
                return True
        else:
            results.append(
                self.get_validation_result(
                    classification, False, 1, allele,
                    self.human_description(t, s),
                    self.concise_description(t, s), errors,
                    gene_tokens, mane_transcript
                )
            )
            return False

    def add_mane_data(self, hgvs_expr, mane, mane_data, s) -> None:
        """Add mane transcript information to mane_data.

        :param str hgvs_expr: HGVS expression of MANE Transcript
        :param dict mane: MANE Transcript information
        :param dict mane_data: MANE Transcript data found for given query
        :param Classification s: Classification token
        """
        key = '_'.join(mane['status'].lower().split())
        if hgvs_expr in mane_data[key].keys():
            mane_data[key][hgvs_expr]['count'] += 1
        else:
            mane_data[key][hgvs_expr] = {
                'classification_token': s,
                'accession': mane['refseq'],
                'count': 1
            }

    def add_mane_to_validation_results(self, mane_data, valid_alleles,
                                       results, classification, gene_tokens):
        """Add MANE Transcript data to list of validation results.

        :param dict mane_data: MANE Transcript data found for given query
        :param list valid_alleles: A list containing current valid alleles
        :param list results: A list of validation results
        :param Classification classification: The classification for tokens
        :param list gene_tokens: List of GeneMatchTokens
        """
        hgvs_exprs = mane_data.keys()
        for key in ['mane_select', 'mane_plus_clinical',
                    'longest_compatible_remaining', 'grch38']:
            highest_count = 0
            mane_result = None
            mane_allele = None
            mane_transcript = None
            if key not in hgvs_exprs:
                continue
            for hgvs_expr in mane_data[key].keys():
                data = mane_data[key][hgvs_expr]
                tmp_allele = None

                if '=' in hgvs_expr:
                    mane_ac = hgvs_expr.split(':')[0]
                    pos = int(re.findall(r'\d+', hgvs_expr.split(':')[-1])[-1])

                    try:
                        sequence_id = self.dp.translate_sequence_identifier(
                            mane_ac, 'ga4gh'
                        )[0]
                    except KeyError:
                        logger.warning(f"GA4GH Data Proxy unable to translate"
                                       f"sequence identifier {mane_ac}")
                    else:
                        new_nuc = self.seqrepo_access.sequence_at_position(
                            mane_ac, pos
                        )
                        if new_nuc:
                            seq_location = models.SequenceLocation(
                                sequence_id=sequence_id,
                                interval=models.SimpleInterval(
                                    start=pos - 1,
                                    end=pos
                                )
                            )

                            state = models.SequenceState(sequence=new_nuc)
                            allele = models.Allele(location=seq_location,
                                                   state=state)
                            allele['_id'] = ga4gh_identify(allele)
                            allele = allele.as_dict()
                            allele_seq_id = allele['location']['sequence_id']
                            if allele_seq_id.startswith('ga4gh:GS.'):
                                allele['location']['sequence_id'] = \
                                    allele_seq_id.replace('ga4gh:GS.',
                                                          'ga4gh:SQ.')

                            if allele:
                                len_of_seq = self.seqrepo_access.len_of_sequence(mane_ac)  # noqa: E501

                                if len_of_seq >= pos - 1:
                                    tmp_allele = allele
                else:
                    tmp_allele = self.get_allele_from_hgvs(hgvs_expr, [])
                if data['count'] > highest_count and tmp_allele:
                    highest_count = data['count']
                    mane_result = data
                    mane_allele = tmp_allele
                    mane_transcript = hgvs_expr

            if mane_allele:
                self.add_validation_result(
                    mane_allele, valid_alleles, results, classification,
                    mane_result['classification_token'],
                    mane_result['accession'], gene_tokens, [],
                    mane_transcript=mane_transcript
                )
                return
