"""Module for Validation."""
from typing import List, Tuple, Optional
from abc import ABC, abstractmethod
from variation.schemas.classification_response_schema import Classification, \
    ClassificationType
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.validation_response_schema import ValidationResult, \
    LookupType
from variation.tokenizers import GeneSymbol
from variation.data_sources import SeqRepoAccess, TranscriptMappings
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
import hgvs.parser
import requests
from requests.exceptions import HTTPError
import logging
from variation.validators.genomic_base import GenomicBase

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class Validator(ABC):
    """The validator class."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol) -> None:
        """Initialize the DelIns validator.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings to/from gene symbols
        :param GeneSymbol gene_symbol: GeneSymbol tokenizer
        """
        self.transcript_mappings = transcript_mappings
        self.seqrepo_access = seqrepo_access
        self._gene_matcher = gene_symbol
        self.dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        self.tlr = Translator(data_proxy=self.dp)
        self.hgvs_parser = hgvs.parser.Parser()
        self.genomic_base = GenomicBase(self.dp)

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
                                  classification, results, gene_tokens) \
            -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of classification Tokens
        :param list transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param list results: Stores validation result objects
        :param list gene_tokens: List of GeneMatchTokens for a classification
        """
        raise NotImplementedError

    def validate(self, classification: Classification) \
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

        self.get_valid_invalid_results(classification_tokens, transcripts,
                                       classification, results, gene_tokens)
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
                              mane_transcript=None) -> None:
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
        else:
            results.append(
                self.get_validation_result(
                    classification, False, 1, allele,
                    self.human_description(t, s),
                    self.concise_description(t, s), errors,
                    gene_tokens, mane_transcript
                )
            )

    def add_mane_transcript(self, classification, results, gene_tokens,
                            mane_transcripts_dict) -> None:
        """Add MANE transcript validation result objects to a list of results.

        :param Classification classification: A classification for a list of
            tokens
        :param list results: A list to store validation result objects
        :param list gene_tokens: List of GeneMatchTokens
        :param dict mane_transcripts_dict: Possible MANE select transcripts
            with classification and transcript
        """
        mane_transcripts = list()
        found_mane_transcripts = list()
        replace_old_keys = list()
        for hgvs_expr in mane_transcripts_dict.keys():
            self.add_valid_mane_transcripts(hgvs_expr, mane_transcripts_dict,
                                            replace_old_keys, gene_tokens,
                                            found_mane_transcripts,
                                            mane_transcripts)

        for tup in replace_old_keys:
            mane_transcripts_dict[tup[1]] = mane_transcripts_dict[tup[0]]
            del mane_transcripts_dict[tup[0]]

        errors = list()

        if len(mane_transcripts) == 0:
            logger.warning("No MANE Select transcript found for "
                           f"{mane_transcripts_dict.keys()}")
            return

        for el in mane_transcripts:
            (hgvs_expr, ensembl_nuc, refseq_nuc,
             ensembl_protein, refseq_protein) = el

            keys = mane_transcripts_dict[hgvs_expr].keys()

            if 'protein' in keys:
                allele = self.get_allele_from_hgvs(refseq_protein, errors)
                if mane_transcripts_dict[hgvs_expr]['protein']:
                    transcript = ensembl_protein
                else:
                    transcript = refseq_protein
            else:
                allele = self.get_allele_from_hgvs(refseq_nuc, errors)
                if mane_transcripts_dict[hgvs_expr]['nucleotide']:
                    transcript = ensembl_nuc
                else:
                    transcript = refseq_nuc

            self.add_validation_result(
                allele, [], results, classification,
                mane_transcripts_dict[hgvs_expr]['classification_token'],
                mane_transcripts_dict[hgvs_expr]['transcript_token'],
                gene_tokens, errors, transcript
            )

    def add_valid_mane_transcripts(self, hgvs_expr, mane_transcripts_dict,
                                   replace_old_keys, gene_tokens,
                                   found_mane_transcripts, mane_transcripts)\
            -> None:
        """Add valid mane transcripts to a list.

        :param str hgvs_expr: The HGVS expression to find a mane transcript for
        :param dict mane_transcripts_dict: Possible MANE select transcripts
            with classification and transcript
        :param list replace_old_keys: Contains tuples of old and new keys to
            replace in mane_transcripts_dict
        :param list gene_tokens: A list of GeneMatchTokens
        :param list found_mane_transcripts: Contains tuple of ensembl and
            refseq hgvs expressions that have already been found
        :param list mane_transcripts: The list of mane transcripts for
            a given hgvs_expr
        """
        if '(' in hgvs_expr and ')' in hgvs_expr:
            # ClinGen Allele Registry doesn't like () in the query
            if mane_transcripts_dict[hgvs_expr][
                'classification_token'].token_type == \
                    'PolypeptideTruncation':  # TODO: Change
                old = hgvs_expr
                hgvs_expr = hgvs_expr.replace('(', '')
                hgvs_expr = hgvs_expr.replace(')', '')
                old_and_new = old, hgvs_expr
                if old_and_new not in replace_old_keys:
                    replace_old_keys.append(old_and_new)
        mane_transcript_tuple = self.get_mane_transcript(hgvs_expr,
                                                         gene_tokens)
        if mane_transcript_tuple:
            ensembl_nuc, refseq_nuc, ensembl_protein, refseq_protein =\
                mane_transcript_tuple

            if ensembl_nuc and refseq_nuc:
                if (ensembl_nuc, refseq_nuc) not in found_mane_transcripts:
                    found_mane_transcripts.append(
                        (ensembl_nuc, refseq_nuc,
                         ensembl_protein, refseq_protein)
                    )
                    mane_transcripts.append(
                        (hgvs_expr, ensembl_nuc, refseq_nuc,
                         ensembl_protein, refseq_protein)
                    )

    def get_mane_transcript(self, hgvs_expr, gene_tokens)\
            -> Optional[Tuple[str, str, str, str]]:
        """Return MANE Select Transcript from ClinGene Allele Registry API.

        :param str hgvs_expr: The HGVS expression to query
        :param list gene_tokens: List of GeneMatchTokens
        :return: Tuple[nucleotide Ensembl HGVS expr, nucleotide RefSeq hgvs
            expr, protein Ensembl hgvs expr, protein RefSeq hgvs expr]
        """
        request = requests.get(
            f"https://reg.genome.network/allele?hgvs={hgvs_expr}")
        if request.status_code == 200:
            resp = request.json()
            if 'transcriptAlleles' in resp.keys():
                mane_transcript = None
                gene_symbol = None
                for t in resp['transcriptAlleles']:
                    if 'geneSymbol' in t.keys() and not gene_symbol:
                        gene_symbol = t['geneSymbol']
                    if 'MANE' in t.keys() and not mane_transcript:
                        mane_transcript = t['MANE']
                    if gene_symbol and mane_transcript:
                        break
                if mane_transcript:
                    if mane_transcript['maneStatus'] == 'MANE Select':
                        if gene_symbol:
                            gene = self._gene_matcher.match(gene_symbol)
                            if gene not in gene_tokens:
                                gene_tokens.append(gene)
                        return (
                            mane_transcript['nucleotide']['Ensembl']['hgvs'],
                            mane_transcript['nucleotide']['RefSeq']['hgvs'],
                            mane_transcript['protein']['Ensembl']['hgvs'],
                            mane_transcript['protein']['RefSeq']['hgvs']
                        )
            else:
                if 'aminoAcidAlleles' in resp.keys() and len(resp['aminoAcidAlleles']) > 0:  # noqa: E501
                    amino_acid_allele = resp['aminoAcidAlleles'][0]
                    if 'hgvsMatchingTranscriptVariant' in amino_acid_allele.keys():  # noqa: E501
                        if len(amino_acid_allele['hgvsMatchingTranscriptVariant']) > 0:  # noqa: E501
                            for t in amino_acid_allele['hgvsMatchingTranscriptVariant']:  # noqa: E501
                                if '[' not in t:
                                    # Temp condition since variation norm
                                    # cant handle multiple possible variations
                                    return self.get_mane_transcript(
                                        t, gene_tokens
                                    )
        return None
