"""Module for Validation."""
from typing import List, Optional, Dict
from abc import ABC, abstractmethod
from variation.schemas.classification_response_schema import Classification, \
    ClassificationType
from variation.schemas.token_response_schema import GeneMatchToken
import variation.schemas.token_response_schema as token_schema
from variation.schemas.validation_response_schema import ValidationResult, \
    LookupType
from variation.tokenizers import GeneSymbol
from variation.data_sources import SeqRepoAccess, TranscriptMappings
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
import hgvs.parser
import logging
from ga4gh.vrs import models, normalize
from ga4gh.core import ga4gh_identify
from variation.validators.genomic_base import GenomicBase
from variation.data_sources import UTA
from bioutils.accessions import coerce_namespace

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
                                  normalize_endpoint, mane_data_found) -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of classification Tokens
        :param list transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param list results: Stores validation result objects
        :param list gene_tokens: List of GeneMatchTokens for a classification
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param dict mane_data_found: MANE Transcript information found
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

        mane_data_found = {
            'mane_select': dict(),
            'mane_plus_clinical': dict(),
            'longest_compatible_remaining': dict(),
            'GRCh38': dict()
        }

        self.get_valid_invalid_results(
            classification_tokens, transcripts, classification,
            results, gene_tokens, normalize_endpoint, mane_data_found
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

    def get_accession(self, t, classification) -> str:
        """Return accession for a classification

        :param str t: Accession
        :param Token classification: Classification for token
        :return: Accession
        """
        if 'HGVS' in classification.matching_tokens or \
                'ReferenceSequence' in classification.matching_tokens:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, token_schema.Token) and t.token_type
                          in ['HGVS', 'ReferenceSequence']][0]
            hgvs_expr = hgvs_token.input_string
            t = hgvs_expr.split(':')[0]
        return t

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

    def to_vrs_allele(self, ac, start, end, coordinate, alt_type, errors,
                      cds_start=None, alt=None) -> Optional[Dict]:
        """Translate accession and position to VRS Allele Object.

        :param str ac: Accession
        :param int start: Start position change
        :param int end: End position change
        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param str alt_type: Type of alteration
        :param list errors: List of errors
        :param int cds_start: Coding start site
        :param str alt: Alteration
        :return: VRS Allele Object
        """
        sequence_id = coerce_namespace(ac)
        try:
            start = int(start)
            if end is None:
                end = start
            end = int(end)
        except (ValueError, TypeError):
            errors.append("Start/End must be valid ints")
            return None

        if coordinate == 'c':
            if cds_start is None:
                cds_start = self.uta.get_cds_start_end(ac)
                if cds_start is None:
                    errors.append(f"Unable to get CDS start for {ac}")
                    return None
            try:
                cds_start = int(cds_start)
            except ValueError:
                errors.append(f"CDS start {cds_start} is not a valid int")
                return None
            start += cds_start
            end += cds_start

        ival_start = start
        ival_end = end

        if alt_type == 'insertion':
            state = alt
            ival_end = ival_start
        elif alt_type in ['substitution', 'deletion', 'delins',
                          'silent_mutation', 'nonsense']:
            if alt_type == 'silent_mutation':
                state = self.seqrepo_access.sequence_at_position(
                    ac, ival_start
                )
            else:
                state = alt or ''
            ival_start -= 1
        else:
            errors.append(f"alt_type not supported: {alt_type}")
            return None

        interval = models.SimpleInterval(start=ival_start, end=ival_end)
        location = models.Location(sequence_id=sequence_id, interval=interval)
        sstate = models.SequenceState(sequence=state)
        allele = models.Allele(location=location, state=sstate)

        try:
            allele = normalize(allele, self.dp)
        except KeyError:
            errors.append(f"Unable to normalize allele: {allele.as_dict()}")
            return None

        if not allele:
            errors.append(f"Unable to find allele for accession, {ac}, "
                          f"and position ({start}, {end})")
            return None

        if coordinate == 'c':
            allele.location.interval.start._value = \
                allele.location.interval.start._value - cds_start
            allele.location.interval.end._value = \
                allele.location.interval.end._value - cds_start

        seq_id = self.dp.translate_sequence_identifier(
            allele.location.sequence_id._value, "ga4gh")[0]
        allele.location.sequence_id = seq_id
        allele._id = ga4gh_identify(allele)
        return allele.as_dict()

    def add_mane_data(self, mane, mane_data, coordinate, alt_type, s,
                      gene_tokens, alt=None) -> None:
        """Add mane transcript information to mane_data.

        :param dict mane: MANE Transcript information
        :param dict mane_data: MANE Transcript data found for given query
        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param str alt_type: Type of alteration
        :param Token s: Classification token
        :param int cds_start: Coding start site
        :param str alt: Alteration
        """
        if not mane:
            return

        if coordinate == 'g':
            if mane['status'] != 'GRCh38':
                s.molecule_context = 'transcript'
                s.reference_sequence = 'c'

                if isinstance(s, token_schema.GenomicSubstitutionToken):
                    # TODO: Only if on different strands
                    ref_rev = s.ref_nucleotide[::-1]
                    alt_rev = s.new_nucleotide[::-1]

                    complements = {
                        'A': 'T',
                        'T': 'A',
                        'C': 'G',
                        'G': 'C'
                    }

                    s.ref_nucleotide = ''
                    s.new_nucleotide = ''
                    for nt in ref_rev:
                        s.ref_nucleotide += complements[nt]
                    for nt in alt_rev:
                        s.new_nucleotide += complements[nt]

        mane_allele = self.to_vrs_allele(
            mane['refseq'], mane['pos'][0], mane['pos'][1],
            coordinate, alt_type, [],
            cds_start=mane.get('coding_start_site', None), alt=alt
        )

        if not mane_allele:
            return

        if coordinate == 'g':
            if not gene_tokens and mane['gene']:
                gene_tokens.append(
                    self._gene_matcher.match(mane['gene'])
                )

        allele_id = mane_allele['_id']

        key = '_'.join(mane['status'].lower().split())
        if allele_id in mane_data[key].keys():
            mane_data[key][allele_id]['count'] += 1
        else:
            mane_data[key][allele_id] = {
                'classification_token': s,
                'accession': mane['refseq'],
                'count': 1,
                'allele': mane_allele,
                'label': mane['refseq']  # TODO: Use VRS to translate
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
        mane_data_keys = mane_data.keys()

        for key in ['mane_select', 'mane_plus_clinical',
                    'longest_compatible_remaining', 'grch38']:
            highest_count = 0
            mane_result = None
            mane_allele = None
            mane_transcript = None
            if key not in mane_data_keys:
                continue
            for mane_allele_id in mane_data[key].keys():
                data = mane_data[key][mane_allele_id]

                if data['count'] > highest_count:
                    highest_count = data['count']
                    mane_result = data
                    mane_allele = data['allele']
                    mane_transcript = data['label']

            if mane_allele:
                self.add_validation_result(
                    mane_allele, valid_alleles, results, classification,
                    mane_result['classification_token'],
                    mane_result['accession'], gene_tokens, [],
                    mane_transcript=mane_transcript
                )
                return
