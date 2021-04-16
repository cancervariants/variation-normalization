"""The module for Single Nucleotide Variant Validation."""
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

        transcripts = self.transcript_mappings.genomic_transcripts(
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

        state = models.SequenceState(sequence=s.new_nucleotide)
        allele = models.Allele(location=seq_location,
                               state=state)
        allele['_id'] = ga4gh_identify(allele)
        return allele.as_dict()

    def get_hgvs_expr(self, classification, t) -> tuple:
        """Return HGVS expression and whether or not it's an Ensembl transcript

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retreived from transcript mapping
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
                valid = True
                errors = list()
                ref_nuc = \
                    self.seq_repo_access.protein_at_position(t, s.position)

                if 'HGVS' in classification.matching_tokens and \
                        not t.startswith('ENST'):
                    hgvs_expr, is_ensembl_transcript = \
                        self.get_hgvs_expr(classification, t)

                    # MANE Select Transcript for HGVS expressions
                    mane_transcripts_dict[hgvs_expr] = {
                        'classification_token': s,
                        'transcript_token': t,
                        'is_ensembl_transcript': is_ensembl_transcript
                    }
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
                    try:
                        sequence_id = \
                            self.dp.translate_sequence_identifier(t,
                                                                  'ga4gh')[0]
                    except KeyError:
                        allele = None
                        error = "GA4GH Data Proxy unable to translate" \
                                f" sequence identifier: {t}"
                        valid = False
                        logger.warning(error)
                        errors.append(error)
                    else:
                        allele = self.get_vrs_allele(sequence_id, s)

                        # MANE Select Transcript for Gene Name + Variation
                        # (ex: BRAF V600E)
                        refseq = ([a for a in self.seq_repo_access.aliases(t)
                                   if a.startswith('refseq:NM_')] or [None])[0]
                        if refseq:
                            gene_token = [t for t in classification.all_tokens
                                          if t.token_type == 'GeneSymbol']
                            if gene_token:
                                is_ensembl_transcript = True
                            else:
                                is_ensembl_transcript = False
                            if 'CodingDNASubstitution' in\
                                    classification.matching_tokens:
                                hgvs_expr = f"{refseq.split('refseq:')[-1]}:" \
                                            f"c.{s.position}" \
                                            f"{s.ref_nucleotide}" \
                                            f">{s.new_nucleotide}"
                                if hgvs_expr not in \
                                        mane_transcripts_dict.keys():
                                    mane_transcripts_dict[hgvs_expr] = {
                                        'classification_token': s,
                                        'transcript_token': t,
                                        'is_ensembl_transcript':
                                            is_ensembl_transcript
                                    }

                if not errors:
                    if ref_nuc != s.ref_nucleotide:
                        errors.append(f'Needed to find {s.ref_nucleotide} at'
                                      f' position {s.position} on {t}'
                                      f' but found {ref_nuc}')
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

        # Now add Mane transcripts to results
        self.add_mane_transcript(classification, results, gene_tokens,
                                 mane_transcripts_dict)

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
            mane_transcript_tuple = self.get_mane_transcript(hgvs_expr)
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

    def get_mane_transcript(self, hgvs_expr):
        """Return MANE Select Transcript from ClinGene Allele Registry.

        :param str hgvs_expr: The HGVS expression to query
        :param bool is_ensembl_transcript: Whether or not the transcript
            is an Ensembl Transcript
        :return: The HGVS RefSeq MANE Select Transcript represented as a string
        """
        request = requests.get(
            f"https://reg.genome.network/allele?hgvs={hgvs_expr}")
        if request.status_code == 200:
            resp = request.json()
            if 'transcriptAlleles' in resp.keys():
                mane_transcript = None
                for t in resp['transcriptAlleles']:
                    if 'MANE' in t.keys():
                        mane_transcript = t['MANE']
                        break
                if mane_transcript:
                    if mane_transcript['maneStatus'] == 'MANE Select':
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
                ([t.token for t in classification.all_tokens if
                 t.token_type in ['HGVS', 'ReferenceSequence']] or [None])[0]

            if not refseq:
                return []

            if ':' in refseq:
                refseq = refseq.split(':')[0]

            res = self.seq_repo_access.aliases(refseq)
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
                                self._gene_matcher.match(gene_symbol))
                    else:
                        logger.warning(f"No gene symbol found for rna "
                                       f"{alias} in transcript_mappings.tsv")
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
        return f'{transcript} {token.ref_nucleotide}' \
               f'{token.position}{token.new_nucleotide}'

    @abstractmethod
    def human_description(self, transcript, token) -> str:
        """Return a human description of the identified variant."""
        raise NotImplementedError