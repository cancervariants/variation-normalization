"""Module for Validation."""
import copy
from typing import List, Optional, Dict, Tuple, Union
from abc import ABC, abstractmethod
from variation.schemas.classification_response_schema import Classification, \
    ClassificationType
from variation.schemas.token_response_schema import GeneMatchToken, Token, \
    GenomicSubstitutionToken
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
from gene.query import QueryHandler as GeneQueryHandler
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class Validator(ABC):
    """The validator class."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTA, dp: SeqRepoDataProxy, tlr: Translator,
                 gene_normalizer: GeneQueryHandler) -> None:
        """Initialize the DelIns validator.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: Gene symbol tokenizer
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param UTA uta: Access to UTA queries
        :param Translator tlr: Translator class
        :param GeneQueryHandler gene_normalizer: Access to gene-normalizer
        """
        self.transcript_mappings = transcript_mappings
        self.seqrepo_access = seqrepo_access
        self._gene_matcher = gene_symbol
        self.dp = dp
        self.tlr = tlr
        self.hgvs_parser = hgvs.parser.Parser()
        self.uta = uta
        self.genomic_base = GenomicBase(self.dp, self.uta)
        self.mane_transcript = mane_transcript
        self.gene_normalizer = gene_normalizer

    @abstractmethod
    def is_token_instance(self, t: Token) -> bool:
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
    def human_description(self, transcript: str, token: Token) -> str:
        """Return a human description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: Human description of the variation change
        """
        raise NotImplementedError

    @abstractmethod
    def concise_description(self, transcript: str, token: Token) -> str:
        """Return a HGVS description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression
        """
        raise NotImplementedError

    @abstractmethod
    def get_gene_tokens(
            self, classification: Classification) -> List[GeneMatchToken]:
        """Return a list of gene tokens for a classification.

        :param Classification classification: Classification for a list of
            tokens
        :return: A list of gene tokens for the classification
        """
        raise NotImplementedError

    @abstractmethod
    def get_transcripts(self, gene_tokens: List,
                        classification: Classification,
                        errors: List) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        raise NotImplementedError

    @abstractmethod
    def validates_classification_type(
            self, classification_type: ClassificationType) -> bool:
        """Check that classification type can be validated by validator.

        :param ClassificationType classification_type: The type of variation
        :return: `True` if classification_type matches validator's
            classification type. `False` otherwise.
        """
        raise NotImplementedError

    @abstractmethod
    def get_valid_invalid_results(
            self, classification_tokens: List, transcripts: List,
            classification: Classification, results: List, gene_tokens: List,
            normalize_endpoint: bool, mane_data_found: Dict,
            is_identifier: bool, hgvs_dup_del_mode: HGVSDupDelModeEnum)\
            -> None:
        """Add validation result objects to a list of results.

        :param List classification_tokens: A list of classification Tokens
        :param List transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param List results: Stores validation result objects
        :param List gene_tokens: List of GeneMatchTokens for a classification
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param Dict mane_data_found: MANE Transcript information found
        :param bool is_identifier: `True` if identifier is given for exact
            location. `False` otherwise.
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `cnv`,
            `repeated_seq_expr`, `literal_seq_expr`.
            This parameter determines how to represent HGVS dup/del expressions
            as VRS objects.
        """
        raise NotImplementedError

    def validate(
            self, classification: Classification, normalize_endpoint: bool,
            hgvs_dup_del_mode: HGVSDupDelModeEnum = HGVSDupDelModeEnum.DEFAULT
    ) -> List[ValidationResult]:
        """Return validation result for a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `cnv`,
            `repeated_seq_expr`, `literal_seq_expr`.
            This parameter determines how to represent HGVS dup/del expressions
            as VRS objects.
        :return: List of ValidationResult's containing valid and invalid
            results
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
            return [
                self.get_validation_result(
                    classification, None, False, 0, {}, '', '',
                    errors, gene_tokens
                )
            ]

        mane_data_found = {
            'mane_select': dict(),
            'mane_plus_clinical': dict(),
            'longest_compatible_remaining': dict(),
            'grch38': dict()
        }

        # If is_identifier, should only run once
        if 'HGVS' in classification.matching_tokens:
            is_identifier = True
        else:
            is_identifier = False

        self.get_valid_invalid_results(
            classification_tokens, transcripts, classification,
            results, gene_tokens, normalize_endpoint, mane_data_found,
            is_identifier, hgvs_dup_del_mode
        )
        return results

    @staticmethod
    def get_validation_result(
            classification: Classification, classification_token: Token,
            is_valid: bool, confidence_score: int, variation: Dict,
            human_description: str, concise_description: str, errors: List,
            gene_tokens: List, identifier: str = None,
            is_mane_transcript: bool = False) -> ValidationResult:
        """Return a validation result object.

        :param Classification classification: The classification for tokens
        :param Token classification_token: Classification token
        :param bool is_valid: Whether or not the classification is valid
        :param int confidence_score: The classification confidence score
        :param Dict variation: A VRS Variation object
        :param str human_description: A human description describing the
            variation
        :param str concise_description: HGVS expression for variation
        :param List errors: A list of errors for the classification
        :param List gene_tokens: List of GeneMatchTokens
        :param str identifier: Identifier for variation
        :param bool is_mane_transcript: `True` if result is MANE transcript.
            `False` otherwise.
        :return: A validation result
        """
        return ValidationResult(
            classification=classification,
            classification_token=classification_token,
            is_valid=is_valid,
            confidence_score=confidence_score,
            variation=variation,
            human_description=human_description,
            concise_description=concise_description,
            errors=errors,
            gene_tokens=gene_tokens,
            is_mane_transcript=is_mane_transcript,
            identifier=identifier
        )

    def get_protein_transcripts(self, gene_tokens: List,
                                errors: List) -> Optional[List[str]]:
        """Get transcripts for variations with protein reference sequence.

        :param List gene_tokens: List of gene tokens for a classification
        :param List errors: List of errors
        :return: List of possible transcript accessions for the variation
        """
        transcripts = self.transcript_mappings.protein_transcripts(
            gene_tokens[0].token, LookupType.GENE_SYMBOL)
        if not transcripts:
            errors.append(f'No transcripts found for gene symbol '
                          f'{gene_tokens[0].token}')
        return transcripts

    def get_coding_dna_transcripts(self, gene_tokens: List,
                                   errors: List) -> Optional[List[str]]:
        """Get transcripts for variations with coding DNA reference sequence.

        :param List gene_tokens: List of gene tokens for a classification
        :param List errors: List of errors
        :return: List of possible transcript accessions for the variation
        """
        transcripts = self.transcript_mappings.coding_dna_transcripts(
            gene_tokens[0].token, LookupType.GENE_SYMBOL)
        if not transcripts:
            errors.append(f'No transcripts found for gene symbol '
                          f'{gene_tokens[0].token}')
        return transcripts

    def get_genomic_transcripts(self, classification: Classification,
                                errors: List) -> Optional[List[str]]:
        """Get NC accessions for variations with genomic reference sequence.

        :param Classification classification: Classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of possible NC accessions for the variation
        """
        nc_accessions = self.genomic_base.get_nc_accessions(classification)
        if not nc_accessions:
            errors.append('Could not find NC_ accession for '
                          f'{self.variation_name()}')
        return nc_accessions

    def get_classification_tokens(
            self, classification: Classification
    ) -> List[Optional[Classification]]:
        """Get classification tokens for a given instance.

        :param Classification classification: A classification for a list of
            tokens
        :return: A list of classification tokens
        """
        return [t for t in classification.all_tokens
                if self.is_token_instance(t)]

    def get_gene_symbol_tokens(
            self, classification: Classification
    ) -> List[Optional[GeneMatchToken]]:
        """Return tokens with GeneSymbol token type from a classification.

        :param Classification classification: Classification of input string
        :return: List of Gene Match Tokens
        """
        return [t for t in classification.all_tokens
                if t.token_type == 'GeneSymbol']

    def _add_gene_symbol_to_tokens(self, gene_symbol: str, gene_symbols: List,
                                   gene_tokens: List) -> None:
        """Add a gene symbol to list of gene match tokens.

        :param str gene_symbol: Gene symbol
        :param List gene_symbols: List of gene symbols matched
        :param List gene_tokens: List of GeneMatchTokens
        """
        if gene_symbol and gene_symbol not in gene_symbols:
            gene_symbols.append(gene_symbol)
            gene_tokens.append(self._gene_matcher.match(
                gene_symbol))

    def _get_gene_tokens(self, classification: Classification,
                         mappings: List) -> List[Optional[GeneMatchToken]]:
        """Get gene symbol tokens for protein or transcript reference
        sequences.

        :param Classification classification: Classification for a list of
            tokens
        :param List mappings: List of transcript mapping methods for
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

    def get_protein_gene_symbol_tokens(
            self, classification: Classification
    ) -> List[Optional[GeneMatchToken]]:
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

    def get_coding_dna_gene_symbol_tokens(
            self, classification: Classification
    ) -> List[Optional[GeneMatchToken]]:
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

    def get_accession(self, t: str, classification: Classification) -> str:
        """Return accession for a classification

        :param str t: Accession
        :param Classification classification: Classification for token
        :return: Accession
        """
        if 'HGVS' in classification.matching_tokens or \
                'ReferenceSequence' in classification.matching_tokens:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type
                          in ['HGVS', 'ReferenceSequence']][0]
            hgvs_expr = hgvs_token.input_string
            t = hgvs_expr.split(':')[0]
        return t

    def add_validation_result(
            self, variation: Dict, valid_variations: List, results: List,
            classification: Classification, s: Token, t: str,
            gene_tokens: List, errors: List, identifier: str = None,
            is_mane_transcript: bool = False) -> bool:
        """Add validation result to list of results.

        :param Dict variation: A VRS Variation object
        :param List valid_variations: A list containing current valid
            variations
        :param List results: A list of validation results
        :param Classification classification: The classification for tokens
        :param Token s: The classification token
        :param string t: Transcript
        :param List gene_tokens: List of GeneMatchTokens
        :param List errors: A list of errors for the classification
        :param str identifier: Identifier for variation
        :param bool is_mane_transcript: `True` if result is MANE transcript.
            `False` otherwise.
        """
        if not errors:
            if is_mane_transcript or \
                    (variation and variation not in valid_variations):
                results.append(
                    self.get_validation_result(
                        classification, s, True, 1, variation,
                        self.human_description(t, s),
                        self.concise_description(t, s), [],
                        gene_tokens,
                        identifier=identifier if identifier else t,
                        is_mane_transcript=is_mane_transcript
                    )
                )
                valid_variations.append(variation)
                return True
        else:
            results.append(
                self.get_validation_result(
                    classification, s, False, 1, variation,
                    self.human_description(t, s),
                    self.concise_description(t, s), errors,
                    gene_tokens,
                    identifier=identifier if identifier else t,
                    is_mane_transcript=is_mane_transcript
                )
            )
            return False

    def _get_start_indef_range(self, start: int) -> models.IndefiniteRange:
        """Return indefinite range given start coordinate

        :param int start: Start position (assumes 1-based)
        :return: Indefinite range model
        """
        return models.IndefiniteRange(value=start - 1, comparator="<=")

    def _get_end_indef_range(self, end: int) -> models.IndefiniteRange:
        """Return indefinite range given end coordinate

        :param int end: End position (assumes 1-based)
        :return: Indefinite range model
        """
        return models.IndefiniteRange(value=end, comparator=">=")

    def _get_ival_certain_range(self, start1: int, start2: int, end1: int,
                                end2: int) -> models.SequenceInterval:
        """Return sequence interval

        :param int start1: Start left pos (assumes 1-based)
        :param int start2: Start right pos (assumes 1-based)
        :param int end1: End left pos (assumes 1-based)
        :param int end2: End right pos (assumes 1-based)
        :return: Sequence Interval model
        """
        return models.SequenceInterval(
            start=models.DefiniteRange(min=start1 - 1, max=start2 - 1),
            end=models.DefiniteRange(min=end1 + 1, max=end2 + 1)
        )

    def _get_sequence_loc(
            self, ac: str, interval: models.SequenceInterval
    ) -> models.Location:
        """Return VRS location

        :param str ac: Accession
        :param models.SequenceInterval interval: VRS sequence interval
        :return: VRS Location model
        """
        return models.Location(sequence_id=coerce_namespace(ac),
                               interval=interval)

    def _get_ival_start_end(
            self, coordinate: str, start: int, end: int, cds_start: int,
            errors: List) -> Optional[Tuple[int, int]]:
        """Get ival_start and ival_end coordinates.

        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param int start: Start position change
        :param int end: End position change
        :param int cds_start: Coding start site
        :param List errors: List of errors
        :return: Tuple[ival_start, ival_end]
        """
        try:
            start = int(start)
            if end is None:
                end = start
            end = int(end)
        except (ValueError, TypeError):
            errors.append("Start/End must be valid ints")
            return None

        if coordinate == 'c':
            if cds_start:
                start += cds_start
                end += cds_start
        return start, end

    def to_vrs_allele_ranges(
            self, ac: str, coordinate: str, alt_type: str, errors: List,
            ival: models.SequenceInterval) -> Optional[Dict]:
        """Translate variation ranges to VRS Allele Object.

        :param str ac: Accession
        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param str alt_type: Type of alteration
        :param list errors: List of errors
        :param models.SequenceInterval ival: Sequence Interval
        :return: VRS Allele object
        """
        if coordinate == 'c':
            # TODO: Once we add support for ranges on c. coord
            return None
        if alt_type in ['uncertain_deletion', 'uncertain_duplication',
                        'duplication_range', 'deletion_range']:
            sstate = models.LiteralSequenceExpression(
                sequence=""
            )
        else:
            errors.append("No state")
            return None

        return self._vrs_allele(ac, ival, sstate, alt_type, errors)

    def _vrs_allele(self, ac: str, interval: models.SequenceInterval,
                    sstate: Union[models.LiteralSequenceExpression,
                                  models.DerivedSequenceExpression,
                                  models.RepeatedSequenceExpression],
                    alt_type: str, errors: List) -> Optional[Dict]:
        """Create a VRS Allele object.

        :param str ac: Accession
        :param SequenceInterval interval: Sequence Interval
        :param sstate: State
        :type sstate: models.LiteralSequenceExpression or
            models.DerivedSequenceExpression or
            models.RepeatedSequenceExpression
        :param str alt_type: Type of alteration
        :param List errors: List of errors
        :return: VRS Allele object represented as a Dict
        """
        location = self._get_sequence_loc(ac, interval)
        allele = models.Allele(location=location, state=sstate)

        # Ambiguous regions do not get normalized
        if alt_type not in ["uncertain_deletion", "uncertain_duplication",
                            "duplication_range", "deletion_range"]:
            try:
                allele = normalize(allele, self.dp)
                if alt_type == 'deletion':
                    allele.state.sequence = ''
            except (KeyError, AttributeError) as e:
                errors.append(f"vrs-python unable to normalize allele: {e}")
                return None

        if not allele:
            errors.append("Unable to get allele")
            return None

        seq_id = self.dp.translate_sequence_identifier(
            allele.location.sequence_id._value, "ga4gh")[0]
        allele.location.sequence_id = seq_id
        allele.location._id = ga4gh_identify(allele.location)
        allele._id = ga4gh_identify(allele)
        return allele.as_dict()

    def to_vrs_allele(
            self, ac: str, start: int, end: int, coordinate: str,
            alt_type: str, errors: List, cds_start: int = None,
            alt: str = None) -> Optional[Dict]:
        """Translate accession and position to VRS Allele Object.

        :param str ac: Accession
        :param int start: Start position change
        :param int end: End position change
        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param str alt_type: Type of alteration
        :param List errors: List of errors
        :param int cds_start: Coding start site
        :param str alt: Alteration
        :return: VRS Allele Object
        """
        ival_coords = self._get_ival_start_end(coordinate, start, end,
                                               cds_start, errors)
        if not ival_coords:
            return None
        if ival_coords[0] > ival_coords[1]:
            ival_end, ival_start = ival_coords
        else:
            ival_start, ival_end = ival_coords

        # Right now, this follows HGVS conventions
        # This will change once we support other representations
        if alt_type == 'insertion':
            state = alt
            ival_end = ival_start
        elif alt_type in ['substitution', 'deletion', 'delins',
                          'silent_mutation', 'nonsense']:
            if alt_type == 'silent_mutation':
                state = self.seqrepo_access.get_sequence(
                    ac, ival_start
                )
                if state is None:
                    errors.append(f"Unable to get sequence on {ac} from "
                                  f"{ival_start}")
                    return None
            else:
                state = alt or ''
            ival_start -= 1
        elif alt_type == 'duplication':
            ref = self.seqrepo_access.get_sequence(ac, ival_start,
                                                   ival_end)
            if ref is not None:
                state = ref + ref
            else:
                errors.append(f"Unable to get sequence on {ac} from "
                              f"{ival_start} to {ival_end}")
                return None
            ival_start -= 1
        else:
            errors.append(f"alt_type not supported: {alt_type}")
            return None

        interval = models.SequenceInterval(
            start=models.Number(value=ival_start),
            end=models.Number(value=ival_end))
        sstate = models.LiteralSequenceExpression(sequence=state)
        return self._vrs_allele(ac, interval, sstate, alt_type, errors)

    def _validate_gene_pos(self, gene: str, alt_ac: str, pos1: int, pos2: int,
                           errors: List, pos3: int = None, pos4: int = None,
                           residue_mode: str = "residue") -> None:
        """Validate whether free text genomic query is valid input.
        If invalid input, add error to list of errors

        :param str gene: Queried gene
        :param str alt_ac: Genomic accession
        :param int pos1: Queried genomic position
        :param int pos2: Queried genomic position
        :param int pos3: Queried genomic position
        :param int pos4: Queried genomic position
        :param str residue_mode: Must be either `inter-residue` or `residue`
        :param List errors: List of errors
        """
        gene_start_end = {"start": None, "end": None}
        resp = self.gene_normalizer.search(gene, incl="Ensembl")
        if resp.source_matches:
            ensembl_resp = resp.source_matches[0]
            if ensembl_resp.records[0].locations:
                ensembl_loc = ensembl_resp.records[0].locations[0]
                gene_start_end["start"] = ensembl_loc.interval.start.value
                gene_start_end["end"] = ensembl_loc.interval.end.value - 1

        if gene_start_end["start"] is None and gene_start_end["end"] is None:
            errors.append(f"gene-normalizer unable to find Ensembl location"
                          f"for {gene}")
        else:
            assembly = self.uta.get_chr_assembly(alt_ac)
            if assembly:
                # Not in GRCh38 assembly. Gene normalizer only uses 38, so we
                # need to liftover to GRCh37 coords
                chromosome, assembly = assembly
                for key in gene_start_end.keys():
                    gene_pos = gene_start_end[key]
                    gene_pos_liftover = \
                        self.uta.liftover_to_37.convert_coordinate(chromosome,
                                                                   gene_pos)
                    if gene_pos_liftover is None or len(gene_pos_liftover) == 0:  # noqa: E501
                        errors.append(f"{gene_pos} does not"
                                      f" exist on {chromosome}")
                        return None
                    else:
                        gene_start_end[key] = gene_pos_liftover[0][1]

            gene_start = gene_start_end["start"]
            gene_end = gene_start_end["end"]

            for pos in [pos1, pos2, pos3, pos4]:
                if pos not in ["?", None]:
                    if residue_mode == "residue":
                        pos -= 1
                    if not (gene_start <= pos <= gene_end):
                        errors.append(f"Position {pos} out of index on "
                                      f"{alt_ac} on gene, {gene}")

    def _get_coord_alt(self, coordinate: str, mane: Dict,
                       s_copy: Token) -> Optional[Tuple[str, str]]:
        """Get coordinate and alteration

        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param Dict mane: Mane data
        :param Token s_copy: classification token
        :return: Coordinate, alteration
        """
        if coordinate == 'g' and mane['status'].lower() != 'grch38':
            s_copy.molecule_context = 'transcript'
            s_copy.reference_sequence = 'c'
            coordinate = s_copy.reference_sequence

            if isinstance(s_copy, GenomicSubstitutionToken) and \
                    mane['strand'] == '-':
                ref_rev = s_copy.ref_nucleotide[::-1]
                alt_rev = s_copy.new_nucleotide[::-1]

                complements = {
                    'A': 'T',
                    'T': 'A',
                    'C': 'G',
                    'G': 'C'
                }

                s_copy.ref_nucleotide = ''
                s_copy.new_nucleotide = ''
                for nt in ref_rev:
                    s_copy.ref_nucleotide += complements[nt]
                for nt in alt_rev:
                    s_copy.new_nucleotide += complements[nt]
                alt = s_copy.new_nucleotide
            else:
                alt = None
            return coordinate, alt
        return None

    def add_mane_data(
            self, mane: Dict, mane_data: Dict, coordinate: str, alt_type: str,
            s: Token, alt: str = None, mane_variation: Dict = None) -> None:
        """Add mane transcript information to mane_data.

        :param Dict mane: MANE data
        :param Dict mane_data: All MANE data found for given query
        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param str alt_type: Type of alteration
        :param Token s: Classification token
        :param str alt: Alteration
        :param Dict mane_variation: VRS Variation for mane data
        """
        if not mane:
            return None

        s_copy = copy.deepcopy(s)
        coord_alt = self._get_coord_alt(coordinate, mane, s_copy)
        if coord_alt:
            coordinate = coord_alt[0] if coord_alt[0] else coordinate
            alt = coord_alt[1] if coord_alt[1] else alt

        if mane_variation is None:
            new_allele = self.to_vrs_allele(
                mane['refseq'], mane['pos'][0], mane['pos'][1],
                coordinate, alt_type, [],
                cds_start=mane.get('coding_start_site', None), alt=alt
            )
            variation = new_allele
        else:
            variation = mane_variation

        if not variation:
            return None

        self._add_dict_to_mane_data(mane['refseq'], s_copy, variation,
                                    mane_data, mane['status'])

    def _add_dict_to_mane_data(self, ac: str, s: Token, variation: Dict,
                               mane_data: Dict, status: str) -> None:
        """Add variation data to mane data for normalize endpoint.

        :param str ac: Accession
        :param Token s: Classification token
        :param Dict variation: VRS Variation object
        :param Dict mane_data: MANE Transcript data found for given query
        :param str status: Status for variation (GRCh38, MANE Select,
            MANE Clinical Plus)
        """
        _id = variation['_id']
        key = '_'.join(status.lower().split())

        if _id in mane_data[key].keys():
            mane_data[key][_id]['count'] += 1
        else:
            mane_data[key][_id] = {
                'classification_token': s,
                'accession': ac,
                'count': 1,
                'variation': variation,
                'label': ac  # TODO: Use VRS to translate
            }

    def add_mane_to_validation_results(
            self, mane_data: Dict, valid_alleles: List, results: List,
            classification: Classification, gene_tokens: List) -> None:
        """Add MANE Transcript data to list of validation results.

        :param Dict mane_data: MANE Transcript data found for given query
        :param List valid_alleles: A list containing current valid alleles
        :param List results: A list of validation results
        :param Classification classification: The classification for tokens
        :param List gene_tokens: List of GeneMatchTokens
        """
        mane_data_keys = mane_data.keys()
        for key in ['mane_select', 'mane_plus_clinical', 'grch38',
                    'longest_compatible_remaining']:
            highest_count = 0
            mane_result = None
            mane_allele = None
            identifier = None
            if key not in mane_data_keys:
                continue
            _mane_data_keys = mane_data[key].keys()
            for mane_allele_id in _mane_data_keys:
                data = mane_data[key][mane_allele_id]

                if data['count'] > highest_count:
                    highest_count = data['count']
                    mane_result = data
                    mane_allele = data['variation']
                    identifier = data['accession']

            if mane_allele:
                self.add_validation_result(
                    mane_allele, valid_alleles, results, classification,
                    mane_result['classification_token'],
                    mane_result['accession'], gene_tokens, [],
                    identifier=identifier, is_mane_transcript=True
                )
                return

    def _check_index(self, ac: str, pos: int, errors: List) -> Optional[str]:
        """Check that index actually exists

        :param str ac: Accession
        :param int pos: Position changes
        :param List errors: List of errors
        :return: Reference sequence
        """
        seq = self.seqrepo_access.get_sequence(ac, pos)
        if not seq:
            errors.append(f"Pos {pos} not found on {ac}")
            return None
        return seq

    def _grch38_dict(self, ac: str, pos: Tuple[int, int]) -> Dict:
        """Create dict for normalized concepts

        :param str ac: Acession
        :param Tuple[int, int] pos: Position changes
        :return: GRCh38 data
        """
        return dict(
            gene=None,
            refseq=ac if ac.startswith('NC') else None,
            ensembl=ac if ac.startswith('ENSG') else None,
            pos=pos,
            strand=None,
            status='GRCh38'
        )

    def _is_grch38_assembly(self, t: str) -> bool:
        """Return whether or not accession is GRCh38 assembly.

        :param str t: Accession
        :return: `True` if accession is GRCh38 assembly. `False` otherwise
        """
        return 'GRCh38' in [a for a in self.dp.get_metadata(t)['aliases'] if a.startswith('GRCh')][0]  # noqa: E501
