"""Module for Validation."""
from typing import List
from abc import ABC, abstractmethod
from variant.schemas.classification_response_schema import Classification, \
    ClassificationType
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.validation_response_schema import ValidationResult
from variant.tokenizers import GeneSymbol
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
import hgvs.parser
import requests
from requests.exceptions import HTTPError
import logging

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class Validator(ABC):
    """The validator class."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol) -> None:
        """Initialize the DelIns validator.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: GeneSymbol tokenizer
        """
        self.transcript_mappings = transcript_mappings
        self.seqrepo_access = seqrepo_access
        self._gene_matcher = gene_symbol
        self.dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        self.tlr = Translator(data_proxy=self.dp)
        self.hgvs_parser = hgvs.parser.Parser()

    @abstractmethod
    def validate(self, classification: Classification) \
            -> List[ValidationResult]:
        """Return validation result for a given classification."""
        raise NotImplementedError

    @abstractmethod
    def validates_classification_type(self,
                                      classification_type: ClassificationType)\
            -> bool:
        """Check that classification type matches."""
        raise NotImplementedError

    @abstractmethod
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
        raise NotImplementedError

    @abstractmethod
    def is_token_instance(self, t) -> bool:
        """Check to see if token is instance of a token type."""
        raise NotImplementedError

    @abstractmethod
    def variant_name(self):
        """Return the variant name"""
        raise NotImplementedError

    @abstractmethod
    def human_description(self, transcript, token) -> str:
        """Return a human description of the identified variant."""
        raise NotImplementedError

    @abstractmethod
    def concise_description(self, transcript, token) -> str:
        """Return a description of the identified variant."""
        raise NotImplementedError

    @abstractmethod
    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return a list of gene tokens for a classification."""
        raise NotImplementedError

    def get_classification_tokens(self, classification):
        """Get classification tokens for a given instance.

        :param Classification classification: A classification for a list of
            tokens
        :return: A list of classification tokens
        """
        return [t for t in classification.all_tokens
                if self.is_token_instance(t)]

    def get_gene_symbol_tokens(self, classification) -> List[GeneMatchToken]:
        """Return tokens with token type GeneSymbol from a classification.

        :param Classification classification: Classification of input string
        :return: List of Gene Match Tokens
        """
        return [t for t in classification.all_tokens
                if t.token_type == 'GeneSymbol']

    def get_protein_gene_symbol_tokens(self, classification)\
            -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        gene_tokens = [t for t in classification.all_tokens
                       if t.token_type == 'GeneSymbol']
        if not gene_tokens:
            # Convert refseq to gene symbol
            refseq = \
                ([t.token for t in classification.all_tokens if
                 t.token_type in ['HGVS', 'ReferenceSequence',
                                  'LocusReferenceGenomic']] or [None])[0]

            if not refseq:
                return []

            if ':' in refseq:
                refseq = refseq.split(':')[0]

            gene_symbols = list()

            for mapping in [
                self.transcript_mappings.get_gene_symbol_from_ensembl_protein,
                self.transcript_mappings.get_gene_symbol_from_refeq_protein,
                self.transcript_mappings.get_gene_symbol_from_lrg
            ]:
                gene_symbol = mapping(refseq)
                self._add_gene_symbol_to_tokens(
                    gene_symbol, gene_symbols, gene_tokens
                )
                if gene_tokens:
                    break

        return gene_tokens

    def _add_gene_symbol_to_tokens(self, gene_symbol, gene_symbols,
                                   gene_tokens):
        """Add a gene symbol to list of gene match tokens.

        :param str gene_symbol: Gene symbol
        :param list gene_symbols: List of gene symbols matched
        :param list gene_tokens: List of GeneMatchTokens
        """
        if gene_symbol and gene_symbol not in gene_symbols:
            gene_symbols.append(gene_symbol)
            gene_tokens.append(self._gene_matcher.match(
                gene_symbol))

    def get_coding_dna_gene_symbol_tokens(self, classification):
        """Return gene symbol tokens for coding dna classifications.

        :param Classification classification: Classification of input string
        :return: A list of gene symbol tokens
        """
        gene_tokens = self.get_gene_symbol_tokens(classification)
        if not gene_tokens:
            # Convert refseq to gene symbol
            refseq = \
                ([t.token for t in classification.all_tokens if
                 t.token_type in ['HGVS', 'ReferenceSequence',
                                  'LocusReferenceGenomic']] or [None])[0]
            if not refseq:
                return []

            if ':' in refseq:
                refseq = refseq.split(':')[0]

            gene_symbols = list()
            for mapping in [
                self.transcript_mappings.get_gene_symbol_from_refseq_rna,
                self.transcript_mappings.get_gene_symbol_from_ensembl_transcript,  # noqa: E501
                self.transcript_mappings.get_gene_symbol_from_lrg
            ]:
                gene_symbol = mapping(refseq)
                self._add_gene_symbol_to_tokens(
                    gene_symbol, gene_symbols, gene_tokens
                )
                if gene_tokens:
                    break
        return gene_tokens

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
            variant
        :param str concise_description: The identified variant
        :param list errors: A list of errors for the classification
        :param list gene_tokens: List of GeneMatchTokens
        :param str mane_transcript: MANE transcript
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

    def get_allele_from_hgvs(self, hgvs_expr, errors):
        """Return allele from a given hgvs_expr.

        :param str hgvs_expr: The HGVS string
        :param list errors: List of errors
        :return: Allele as a dictionary
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
        return allele

    def add_validation_result(self, allele, valid_alleles, results,
                              classification, s, t, gene_tokens, errors,
                              mane_transcript=None):
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
                            mane_transcripts_dict):
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
            hgvs_expr = el[0]
            ensembl_transcript = el[1]
            refseq_transcript = el[2]
            refseq_protein = el[3]

            if 'is_ensembl_transcript' in \
                    mane_transcripts_dict[hgvs_expr].keys():
                allele = self.get_allele_from_hgvs(refseq_transcript, errors)
                if mane_transcripts_dict[hgvs_expr]['is_ensembl_transcript']:
                    transcript = ensembl_transcript
                else:
                    transcript = refseq_transcript
            else:
                transcript = refseq_protein
                allele = self.get_allele_from_hgvs(transcript, errors)

            self.add_validation_result(
                allele, [], results, classification,
                mane_transcripts_dict[hgvs_expr]['classification_token'],
                mane_transcripts_dict[hgvs_expr]['transcript_token'],
                gene_tokens, errors, transcript
            )

    def add_valid_mane_transcripts(self, hgvs_expr, mane_transcripts_dict,
                                   replace_old_keys, gene_tokens,
                                   found_mane_transcripts, mane_transcripts):
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
            ensembl, refseq, refseq_protein = mane_transcript_tuple

            if ensembl and refseq:
                if (ensembl, refseq) not in found_mane_transcripts:
                    found_mane_transcripts.append((ensembl, refseq,
                                                   refseq_protein))
                    mane_transcripts.append((hgvs_expr, ensembl, refseq,
                                             refseq_protein))

    def get_mane_transcript(self, hgvs_expr, gene_tokens):
        """Return MANE Select Transcript from ClinGene Allele Registry.

        :param str hgvs_expr: The HGVS expression to query
        :param list gene_tokens: List of GeneMatchTokens
        :return: The HGVS RefSeq MANE Select Transcript represented as a string
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
                            mane_transcript['protein']['RefSeq']['hgvs']
                        )
            else:
                if 'aminoAcidAlleles' in resp.keys() and len(resp['aminoAcidAlleles']) > 0:  # noqa: E501
                    amino_acid_allele = resp['aminoAcidAlleles'][0]
                    if 'hgvsMatchingTranscriptVariant' in amino_acid_allele.keys():  # noqa: E501
                        if len(amino_acid_allele['hgvsMatchingTranscriptVariant']) > 0:  # noqa: E501
                            for t in amino_acid_allele['hgvsMatchingTranscriptVariant']:  # noqa: E501
                                if '[' not in t:
                                    # Temp condition since variant norm
                                    # cannot handle multiple possible variants
                                    return self.get_mane_transcript(
                                        t, gene_tokens
                                    )
        return None
