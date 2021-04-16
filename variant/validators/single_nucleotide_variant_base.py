"""The module for Single Nucleotide Variant Validation."""
from typing import List
from abc import abstractmethod
from requests.exceptions import HTTPError
from .validator import Validator
from variant.schemas.classification_response_schema import Classification
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.validation_response_schema import ValidationResult
from variant.tokenizers import GeneSymbol
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
import hgvs.parser
import requests
import logging

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class SingleNucleotideVariantBase(Validator):
    """The Single Nucleotide Variant Validator Base class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol) \
            -> None:
        """Initialize the validator.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        """
        self.transcript_mappings = transcript_mappings
        self.seq_repo_access = seq_repo_access
        self._gene_matcher = gene_symbol
        self.dp = SeqRepoDataProxy(seq_repo_access.seq_repo_client)
        self.tlr = Translator(data_proxy=self.dp)
        self.hgvs_parser = hgvs.parser.Parser()

    @abstractmethod
    def validate(self, classification: Classification) \
            -> List[ValidationResult]:
        """Validate a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :return: A list of validation results
        """
        raise NotImplementedError

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

        state = models.SequenceState(sequence=s.new_nucleotide)
        allele = models.Allele(location=seq_location,
                               state=state)
        allele['_id'] = ga4gh_identify(allele)
        return allele.as_dict()

    @abstractmethod
    def get_hgvs_expr(self, classification, t):
        """Return HGVS expression and whether or not it's an Ensembl transcript

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        """
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

    def add_mane_transcript(self, classification, results, gene_tokens,
                            mane_transcripts_dict):
        """Add MANE transcript validation result objects to a list of results.

        :param Classification classification: A classification for a list of
            tokens
        :param list results: A list to store validation result objects
        :param list gene_tokens: List of GeneMatchTokens
        :param dict mane_transcripts_dict: Possible MANE select transcripts
            with classification and Ensembl transcript
        """
        mane_transcripts = list()
        found_mane_transcripts = list()
        replace_old_keys = list()
        for hgvs_expr in mane_transcripts_dict.keys():
            if '(' in hgvs_expr and ')' in hgvs_expr:
                # ClinGen Allele Registry doesn't like () in the query
                if mane_transcripts_dict[hgvs_expr][
                    'classification_token'].token_type ==\
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
                ensembl, refseq = mane_transcript_tuple

                if ensembl and refseq:
                    if (ensembl, refseq) not in found_mane_transcripts:
                        found_mane_transcripts.append((ensembl, refseq))
                        mane_transcripts.append((hgvs_expr, ensembl, refseq))

        for tup in replace_old_keys:
            mane_transcripts_dict[tup[1]] = mane_transcripts_dict[tup[0]]
            del mane_transcripts_dict[tup[0]]

        errors = list()
        allele = {}
        seq_location = None

        if len(mane_transcripts) == 0:
            logger.warning("No MANE Select transcript found for "
                           f"{mane_transcripts_dict.keys()}")
            return

        for el in mane_transcripts:
            hgvs_expr = el[0]
            ensembl_transcript = el[1]
            refseq_transcript = el[2]

            try:
                seq_location = self.tlr.translate_from(refseq_transcript,
                                                       'hgvs')
            except HTTPError:
                errors.append(f"{mane_transcripts[0]} is an invalid HGVS "
                              f"expression.")
            except KeyError:
                allele = None
                error = "GA4GH Data Proxy unable to translate  sequence " \
                        f"identifier: {refseq_transcript}"
                logger.warning(error)
                errors.append(error)
            else:
                allele = seq_location.as_dict()

            if mane_transcripts_dict[hgvs_expr]['is_ensembl_transcript']:
                transcript = ensembl_transcript
            else:
                transcript = refseq_transcript

            if not errors:
                results.append(self.get_validation_result(
                    classification, True, 1, allele,
                    self.human_description(
                        mane_transcripts_dict[hgvs_expr][
                            'transcript_token'],
                        mane_transcripts_dict[hgvs_expr][
                            'classification_token']
                    ),
                    self.concise_description(
                        mane_transcripts_dict[hgvs_expr][
                            'transcript_token'],
                        mane_transcripts_dict[hgvs_expr][
                            'classification_token']
                    ), errors, gene_tokens, transcript
                ))
            else:
                results.append(self.get_validation_result(
                    classification, False, 1, allele,
                    self.human_description(
                        mane_transcripts_dict[hgvs_expr]['transcript_token'],
                        mane_transcripts_dict[hgvs_expr][
                            'classification_token']
                    ),
                    self.concise_description(
                        mane_transcripts_dict[hgvs_expr]['transcript_token'],
                        mane_transcripts_dict[hgvs_expr][
                            'classification_token']
                    ), errors, gene_tokens, transcript
                ))

    def get_mane_transcript(self, hgvs_expr, gene_tokens):
        """Return MANE Select Transcript from ClinGene Allele Registry.

        :param str hgvs_expr: The HGVS expression to query
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
                    if 'MANE' in t.keys():
                        mane_transcript = t['MANE']
                        gene_symbol = t['geneSymbol']
                        break
                if mane_transcript:
                    if mane_transcript['maneStatus'] == 'MANE Select':
                        if gene_symbol:
                            gene = self._gene_matcher.match(gene_symbol)
                            if gene not in gene_tokens:
                                gene_tokens.append(gene)
                        return (mane_transcript['nucleotide']['Ensembl']['hgvs'],   # noqa: E501
                                mane_transcript['nucleotide']['RefSeq']['hgvs'])  # noqa: E501
            else:
                if 'aminoAcidAlleles' in resp.keys() and len(resp['aminoAcidAlleles']) > 0:  # noqa: E501
                    amino_acid_allele = resp['aminoAcidAlleles'][0]
                    if 'hgvsMatchingTranscriptVariant' in amino_acid_allele.keys():  # noqa: E501
                        if len(amino_acid_allele['hgvsMatchingTranscriptVariant']) > 0:  # noqa: E501
                            for t in amino_acid_allele['hgvsMatchingTranscriptVariant']:  # noqa: E501
                                if '>' in t and '[' not in t:
                                    # Temp condition since variant norm
                                    # cannot handle multiple possible variants
                                    return self.get_mane_transcript(t)
        return None

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

    @abstractmethod
    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return a list of gene tokens for a classification."""
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
    def validates_classification_type(self, classification_type) \
            -> bool:
        """Return the classification type"""
        raise NotImplementedError

    def concise_description(self, transcript, token) -> str:
        """Return a description of the identified variant."""
        return f'{transcript} {token.ref_nucleotide}' \
               f'{token.position}{token.new_nucleotide}'

    @abstractmethod
    def human_description(self, transcript, token) -> str:
        """Return a human description of the identified variant."""
        raise NotImplementedError
