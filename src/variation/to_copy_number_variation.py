"""Module for to copy number variation translation"""
import datetime
from typing import Dict, List, NamedTuple, Optional, Tuple, Union
from urllib.parse import unquote

from cool_seq_tool.handlers import SeqRepoAccess
from cool_seq_tool.mappers import LiftOver
from cool_seq_tool.schemas import Assembly
from cool_seq_tool.sources import UtaDatabase
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models
from gene.query import QueryHandler as GeneQueryHandler
from gene.schemas import MatchType as GeneMatchType
from pydantic import ValidationError

from variation.classify import Classify
from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.copy_number_schema import (
    AmplificationToCxVarQuery,
    AmplificationToCxVarService,
    Comparator,
    ParsedPosType,
    ParsedToCnVarQuery,
    ParsedToCnVarService,
    ParsedToCxVarQuery,
    ParsedToCxVarService,
)
from variation.schemas.hgvs_to_copy_number_schema import (
    HgvsToCopyNumberChangeService,
    HgvsToCopyNumberCountService,
)
from variation.schemas.normalize_response_schema import (
    HGVSDupDelModeOption,
    ServiceMeta,
)
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.token_response_schema import TokenType
from variation.schemas.validation_response_schema import ValidationResult
from variation.to_vrs import ToVRS
from variation.tokenize import Tokenize
from variation.translate import Translate
from variation.utils import get_priority_sequence_location
from variation.validate import Validate
from variation.version import __version__

VALID_CLASSIFICATION_TYPES = [
    ClassificationType.GENOMIC_DUPLICATION,
    ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS,
    ClassificationType.GENOMIC_DELETION,
    ClassificationType.GENOMIC_DELETION_AMBIGUOUS,
]


class ToCopyNumberError(Exception):
    """Custom exceptions when representing copy number"""


class ParsedAccessionSummary(NamedTuple):
    """Represents accession for parsed endpoints"""

    accession: str
    lifted_over: bool


class ParsedChromosomeSummary(NamedTuple):
    """Represents chromosome and assembly for parsed endpoints"""

    accession: str
    chromosome: str
    lifted_over: bool


class ToCopyNumberVariation(ToVRS):
    """Class for representing copy number variation"""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        tokenizer: Tokenize,
        classifier: Classify,
        validator: Validate,
        translator: Translate,
        gene_normalizer: GeneQueryHandler,
        uta: UtaDatabase,
        liftover: LiftOver
    ) -> None:
        """Initialize theToCopyNumberVariation class

        :param seqrepo_access: Access to SeqRepo
        :param tokenizer: Instance for tokenizing input strings
        :param classifier: Instance for classifying list of ordered tokens
        :param validator: Instance for validating classification
        :param translator: Instance for translating valid results to VRS representations
        :param gene_normalizer: Client for normalizing gene concepts
        :param uta: Access to UTA queries
        :param liftover: Instance to provide mapping between human genome assemblies
        """
        super().__init__(seqrepo_access, tokenizer, classifier, validator, translator)
        self.gene_normalizer = gene_normalizer
        self.uta = uta
        self.liftover = liftover

    async def _get_valid_results(self, q: str) -> Tuple[List[ValidationResult], List]:
        """Get valid results for to copy number variation endpoint

        :param q: Input query string
        :return: Valid results and list of warnings
        """
        valid_results = []
        warnings = []

        # Get tokens for input query
        tokens = self.tokenizer.perform(unquote(q.strip()), warnings)
        if not tokens:
            return valid_results, warnings

        # Get classification for list of tokens
        classification = self.classifier.perform(tokens)
        if not classification:
            warnings.append(f"Unable to find classification for: {q}")
            return valid_results, warnings

        # Ensure that classification is HGVS duplication or deletion
        tmp_classification = None
        if all(
            (
                classification.classification_type in VALID_CLASSIFICATION_TYPES,
                TokenType.HGVS
                in {t.token_type for t in classification.matching_tokens},
            )
        ):
            tmp_classification = classification

        classification = tmp_classification
        if not classification:
            warnings = [f"{q} is not a supported HGVS genomic duplication or deletion"]
            return valid_results, warnings

        # Get validation summary for classification
        validation_summary = await self.validator.perform(classification)
        if validation_summary.valid_results:
            valid_results = validation_summary.valid_results
        else:
            warnings = validation_summary.warnings
            valid_results = []

        return valid_results, warnings

    async def _hgvs_to_cnv_resp(
        self,
        copy_number_type: HGVSDupDelModeOption,
        do_liftover: bool,
        valid_results: Tuple[List[ValidationResult], List[str]],
        warnings: List[str],
        baseline_copies: Optional[int] = None,
        copy_change: Optional[models.CopyChange] = None,
    ) -> Tuple[
        Optional[Union[models.CopyNumberCount, models.CopyNumberChange]], List[str]
    ]:
        """Return copy number variation and warnings response

        :param copy_number_type: The type of copy number variation. Must be either
            `copy_number_count` or `copy_number_change`
        :param hgvs_expr: HGVS expression
        :param do_liftover: Whether or not to liftover to GRCh38 assembly
        :param Valid results and warnings for hgvs_expr
        :param warnings: List of warnings
        :return: CopyNumberVariation and warnings
        """
        variation = None
        if valid_results:
            if copy_number_type == HGVSDupDelModeOption.COPY_NUMBER_CHANGE:
                endpoint_name = Endpoint.HGVS_TO_COPY_NUMBER_CHANGE
            else:
                endpoint_name = Endpoint.HGVS_TO_COPY_NUMBER_COUNT

            translations, warnings = await self.get_translations(
                valid_results,
                warnings,
                hgvs_dup_del_mode=copy_number_type,
                endpoint_name=endpoint_name,
                copy_change=copy_change,
                baseline_copies=baseline_copies,
                do_liftover=do_liftover,
            )
            if translations:
                variation = translations[0].vrs_variation

        if variation:
            if copy_number_type == HGVSDupDelModeOption.COPY_NUMBER_COUNT:
                variation = models.CopyNumberCount(**variation)
            else:
                variation = models.CopyNumberChange(**variation)
        return variation, warnings

    async def hgvs_to_copy_number_count(
        self,
        hgvs_expr: str,
        baseline_copies: int,
        do_liftover: bool = False,
    ) -> HgvsToCopyNumberCountService:
        """Given hgvs, return abolute copy number variation

        :param hgvs_expr: HGVS expression
        :param baseline_copies: Baseline copies number
        :param do_liftover: Whether or not to liftover to GRCh38 assembly
        :return: HgvsToCopyNumberCountService containing Copy Number Count
            Variation and warnings
        """
        valid_results, warnings = await self._get_valid_results(hgvs_expr)
        cn_var, warnings = await self._hgvs_to_cnv_resp(
            HGVSDupDelModeOption.COPY_NUMBER_COUNT,
            do_liftover,
            valid_results,
            warnings,
            baseline_copies=baseline_copies,
        )

        return HgvsToCopyNumberCountService(
            hgvs_expr=hgvs_expr,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
            ),
            copy_number_count=cn_var,
        )

    async def hgvs_to_copy_number_change(
        self,
        hgvs_expr: str,
        copy_change: Optional[models.CopyChange],
        do_liftover: bool = False,
    ) -> HgvsToCopyNumberChangeService:
        """Given hgvs, return copy number change variation

        :param hgvs_expr: HGVS expression
        :param copy_change: The copy change
        :param do_liftover: Whether or not to liftover to GRCh38 assembly
        :return: HgvsToCopyNumberChangeService containing Copy Number Change
            Variation and warnings
        """
        valid_results, warnings = await self._get_valid_results(hgvs_expr)
        cx_var, warnings = await self._hgvs_to_cnv_resp(
            HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
            do_liftover,
            valid_results,
            warnings,
            copy_change=copy_change,
        )

        return HgvsToCopyNumberChangeService(
            hgvs_expr=hgvs_expr,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
            ),
            copy_number_change=cx_var,
        )

    def _get_parsed_ac(
        self, assembly: ClinVarAssembly, chromosome: str, use_grch38: bool = False
    ) -> ParsedAccessionSummary:
        """Get accession for parsed components

        :param assembly: Assembly
        :param chromosome: Chromosome
        :param use_grch38: Whether or not to use GRCh38 assembly
        :raises ToCopyNumberError: If unable to translate assembly and chromosome
            to an accession
        :return: ParsedAccessionSummary containing accession and whether or not it was
            lifted over
        """
        accession = None
        lifted_over = False
        og_assembly = assembly

        if assembly == ClinVarAssembly.HG38:
            assembly = ClinVarAssembly.GRCH38
        elif assembly == ClinVarAssembly.HG19:
            assembly = ClinVarAssembly.GRCH37
        elif assembly == ClinVarAssembly.HG18:
            assembly = ClinVarAssembly.NCBI36

        if use_grch38 and assembly != ClinVarAssembly.GRCH38:
            lifted_over = True
            assembly = ClinVarAssembly.GRCH38

        if assembly != ClinVarAssembly.NCBI36:
            # Variation Normalizer does not support NCBI36 yet
            query = f"{assembly.value}:{chromosome}"
            aliases, error = self.seqrepo_access.translate_identifier(query, "ga4gh")
            if aliases:
                accession = aliases[0]
            else:
                raise ToCopyNumberError(str(error))
        else:
            msg = f"{og_assembly.value} assembly is not currently supported"
            raise ToCopyNumberError(msg)

        return ParsedAccessionSummary(lifted_over=lifted_over, accession=accession)

    def _get_parsed_ac_chr(
        self, accession: str, do_liftover: bool
    ) -> ParsedChromosomeSummary:
        """Get accession and chromosome for parsed components

        :param accession: Genomic accession
        :param do_liftover: Whether or not to liftover to GRCh38 assembly
        :raises ToCopyNumberError: If unable to translate accession
        :return: ParsedChromosomeSummary containing chromosome, accession, and whether
            or not it was lifted over
        """
        chromosome = None
        new_ac = None
        lifted_over = False

        aliases, error = self.seqrepo_access.translate_identifier(accession)
        if error:
            raise ToCopyNumberError(error)

        grch_aliases = [
            a for a in aliases if a.startswith(("GRCh38:chr", "GRCh37:chr"))
        ]

        if grch_aliases:
            grch_record = grch_aliases[0]
            chromosome = grch_record.split(":")[-1]

            if grch_record.startswith("GRCh38") or not do_liftover:
                new_ac = next(a for a in aliases if a.startswith("ga4gh"))
            else:
                grch38_query = grch_record.replace("GRCh37", "GRCh38")
                aliases, error = self.seqrepo_access.translate_identifier(
                    grch38_query, "ga4gh"
                )

                if error:
                    raise ToCopyNumberError(error)

                lifted_over = True
                new_ac = aliases[0]
        else:
            msg = f"Not a supported genomic accession: {accession}"
            raise ToCopyNumberError(msg)

        return ParsedChromosomeSummary(
            accession=new_ac, chromosome=chromosome, lifted_over=lifted_over
        )

    def _validate_ac_pos(self, accession: str, pos: int) -> None:
        """Validate position for parsed components

        :param accession: Genomic accession
        :param pos: Position on accession
        :raises ToCopyNumberError: If position is not valid on accession or
            if accession is not found in seqrepo
        """
        try:
            ref = self.seqrepo_access.sr[accession][pos - 1]
        except ValueError as e:
            msg = f"SeqRepo ValueError: {str(e).replace('start', 'Position')}"
            raise ToCopyNumberError(msg) from e
        except KeyError as e:
            msg = f"Accession not found in SeqRepo: {accession}"
            raise ToCopyNumberError(msg) from e
        else:
            if ref == "":
                msg = f"Position ({pos}) is not valid on {accession}"
                raise ToCopyNumberError(msg) from None

    def _get_vrs_loc_start_or_end(
        self,
        accession: str,
        pos0: int,
        pos_type: ParsedPosType,
        is_start: bool = True,
        pos1: Optional[int] = None,
        comparator: Optional[Comparator] = None,
    ) -> Union[int, models.Range]:
        """Get VRS Sequence Location start and end values

        :param accession: Genomic accession for sequence
        :param pos0: Position (residue coords). If `pos_type` is a definite range,
            this will be the min start position
        :param pos_type: Type of the pos value in VRS Sequence Location
        :param is_start: `True` if position(s) describing VRS start value. `False` if
            position(s) describing VRS end value
        :param pos1: Only set when end is a definite range, this will be the max end
            position
        :param comparator: Must provide when `pos_type` is an Indefinite Range.
            Indicates which direction the range is indefinite. To represent (#_?), set
            to '<='. To represent (?_#), set to '>='.
        :raises ToCopyNumberError: If position is not valid on accession when
            using definite range
        :return: VRS start or end value for sequence location
        """
        if pos_type == ParsedPosType.NUMBER:
            vrs_val = pos0 - 1 if is_start else pos0
        elif pos_type == ParsedPosType.DEFINITE_RANGE:
            self._validate_ac_pos(accession, pos1)
            vrs_val = models.Range(
                [pos0 - 1 if is_start else pos0, pos1 - 1 if is_start else pos1]
            )
        else:
            if comparator == Comparator.LT_OR_EQUAL:
                vrs_val = models.Range([None, pos0 - 1 if is_start else pos0])
            else:
                vrs_val = models.Range([pos0 - 1 if is_start else pos0, None])

        return vrs_val

    def _get_parsed_seq_loc(
        self,
        accession: str,
        chromosome: str,
        start0: int,
        start_pos_type: ParsedPosType,
        end0: int,
        end_pos_type: ParsedPosType,
        start1: Optional[int] = None,
        end1: Optional[int] = None,
        liftover_pos: bool = False,
        start_pos_comparator: Optional[Comparator] = None,
        end_pos_comparator: Optional[Comparator] = None,
    ) -> Tuple[Optional[Dict], Optional[str]]:
        """Get sequence location for parsed components. Accession will be validated.

        :param accession: Genomic accession for sequence
        :param chromosome: Chromosome
        :param start0: Start position (residue coords). If start is a definite range,
            this will be the min start position
        :param start_pos_type: Type of the start value in VRS Sequence Location
        :param end0: End position (residue coords). If end is a definite range, this
            will be the min end position
        :param end_pos_type: Type of the end value in VRS Sequence Location
        :param start1: Only set when start is a definite range, this will be the max
            start position
        :param end1: Only set when end is a definite range, this will be the max end
            position
        :param liftover_pos: Whether or not to liftover positions
        :param start_pos_comparator: Must provide when `start_pos_type` is an Indefinite
            Range. Indicates which direction the range is indefinite. To represent
            (#_?), set to '<='. To represent (?_#), set to '>='.
        :param end_pos_comparator: Must provide when `end_pos_type` is an Indefinite
            Range. Indicates which direction the range is indefinite. To represent
            (#_?), set to '<='. To represent (?_#), set to '>='.
        :raises ToCopyNumberError: If error lifting over positions, translating
            accession, positions not valid on accession,
        :return: Tuple containing VRS sequence location represented as dict (if valid)
            and warning (if invalid)
        """
        seq_loc = None

        # Liftover pos if needed
        if liftover_pos:
            liftover_pos = self._liftover_pos(chromosome, start0, end0, start1, end1)
            start0 = liftover_pos["start0"]
            end0 = liftover_pos["end0"]
            start1 = liftover_pos["start1"]
            end1 = liftover_pos["end1"]

        sequences, error = self.seqrepo_access.translate_identifier(accession, "ga4gh")
        if error:
            raise ToCopyNumberError(error)

        sequence = sequences[0].split("ga4gh:")[-1]

        for pos in [start0, end0]:
            # validate start0 and end0 since they're always required
            self._validate_ac_pos(accession, pos)

        start_vrs = self._get_vrs_loc_start_or_end(
            accession,
            start0,
            start_pos_type,
            is_start=True,
            pos1=start1,
            comparator=start_pos_comparator,
        )

        end_vrs = self._get_vrs_loc_start_or_end(
            accession,
            end0,
            end_pos_type,
            is_start=False,
            pos1=end1,
            comparator=end_pos_comparator,
        )

        seq_loc = models.SequenceLocation(
            sequenceReference=models.SequenceReference(refgetAccession=sequence),
            start=start_vrs,
            end=end_vrs,
        )
        seq_loc.id = ga4gh_identify(seq_loc)

        return seq_loc.model_dump(exclude_none=True) if seq_loc else seq_loc

    def _liftover_pos(
        self,
        chromosome: str,
        start0: int,
        end0: int,
        start1: Optional[int],
        end1: Optional[int],
    ) -> Dict:
        """Liftover GRCh37 positions to GRCh38 positions

        :param chromosome: Chromosome. Must be contain 'chr' prefix, i.e 'chr7'.
        :param start0: Start position (residue coords) GRCh37 assembly. If start is a
            definite range, this will be the min start position
        :param end0: End position (residue coords) GRCh37 assembly. If end is a definite
            range, this will be the min end position
        :param start1: Only set when start is a definite range, this will be the max
            start position. GRCh37 assembly
        :param end1: Only set when end is a definite range, this will be the max end
            position. GRCh37 assembly
        :raises ToCopyNumberError: If unable to liftover position
        :return: Dictionary containing lifted over positions
            ('start0', 'end0', 'start1', 'end1')
        """
        liftover_pos = {"start0": None, "end0": None, "start1": None, "end1": None}

        for k, pos in [
            ("start0", start0),
            ("end0", end0),
            ("start1", start1),
            ("end1", end1),
        ]:
            if pos is not None:
                liftover = self.liftover.get_liftover(
                    chromosome, pos, Assembly.GRCH38
                )
                if not liftover:
                    msg = f"Unable to liftover: {chromosome} with pos {pos}"
                    raise ToCopyNumberError(msg)

                liftover_pos[k] = liftover[1]

        return liftover_pos

    def parsed_to_copy_number(
        self, request_body: Union[ParsedToCnVarQuery, ParsedToCxVarQuery]
    ) -> Union[ParsedToCnVarService, ParsedToCxVarService]:
        """Given parsed genomic components, return Copy Number Count or Copy Number
        Change Variation

        :param request_body: request body
        :return: If `copy_number_type` is Copy Number Count, return ParsedToCnVarService
            containing Copy Number Count variation and list of warnings. Else, return
            ParsedToCxVarService containing Copy Number Change variation and list of
            warnings
        """
        variation = None
        warnings = []

        is_cx = isinstance(request_body, ParsedToCxVarQuery)
        lifted_over = False

        try:
            if not request_body.accession:
                accession_summary = self._get_parsed_ac(
                    request_body.assembly,
                    request_body.chromosome,
                    use_grch38=request_body.do_liftover,
                )
                chromosome = request_body.chromosome
                accession = accession_summary.accession
                lifted_over = accession_summary.lifted_over
            else:
                chr_summary = self._get_parsed_ac_chr(
                    request_body.accession, request_body.do_liftover
                )
                accession = chr_summary.accession
                chromosome = chr_summary.chromosome
                lifted_over = chr_summary.lifted_over

            seq_loc = self._get_parsed_seq_loc(
                accession,
                chromosome,
                request_body.start0,
                request_body.start_pos_type,
                request_body.end0,
                request_body.end_pos_type,
                start1=request_body.start1,
                end1=request_body.end1,
                start_pos_comparator=request_body.start_pos_comparator,
                end_pos_comparator=request_body.end_pos_comparator,
                liftover_pos=request_body.do_liftover and lifted_over,
            )
        except ToCopyNumberError as e:
            warnings.append(str(e))
        else:
            if is_cx:
                variation = models.CopyNumberChange(
                    location=seq_loc, copyChange=request_body.copy_change
                )
                variation.id = ga4gh_identify(variation)
            else:
                if request_body.copies_type == ParsedPosType.NUMBER:
                    copies = request_body.copies0
                elif request_body.copies_type == ParsedPosType.DEFINITE_RANGE:
                    copies = models.Range([request_body.copies0, request_body.copies1])
                else:
                    if request_body.copies_comparator == Comparator.LT_OR_EQUAL:
                        copies = models.Range([None, request_body.copies0])
                    else:
                        copies = models.Range([request_body.copies0, None])
                variation = models.CopyNumberCount(location=seq_loc, copies=copies)
                variation.id = ga4gh_identify(variation)

        service_params = {
            "warnings": warnings,
            "service_meta_": ServiceMeta(
                version=__version__,
                response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
            ),
        }

        if is_cx:
            service_params["copy_number_change"] = variation
        else:
            service_params["copy_number_count"] = variation

        return (
            ParsedToCxVarService(**service_params)
            if is_cx
            else ParsedToCnVarService(**service_params)
        )

    def amplification_to_cx_var(
        self,
        gene: str,
        sequence_id: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> AmplificationToCxVarService:
        """Return Copy Number Change Variation for Amplification query
        Parameter priority:
            1. sequence_id, start, end (must provide ALL)
            2. use the gene-normalizer to get the SequenceLocation

        :param gene: Gene query
        :param sequence_id: Sequence ID for the location. If set, must also provide
            `start` and `end`
        :param start: Start position as residue coordinate for the sequence location. If
            set, must also provide `sequence` and `end`
        :param end: End position as residue coordinate for the sequence location. If
            set, must also provide `sequence` and `start`
        :return: AmplificationToCxVarService containing Copy Number Change and
            list of warnings
        """
        warnings = []
        amplification_label = None
        variation = None
        try:
            og_query = AmplificationToCxVarQuery(
                gene=gene, sequence_id=sequence_id, start=start, end=end
            )
        except ValidationError as e:
            warnings.append(str(e))
            og_query = None
        else:
            # Need to validate the input gene
            gene_norm_resp = self.gene_normalizer.normalize(gene)
            if gene_norm_resp.match_type != GeneMatchType.NO_MATCH:
                vrs_location = None
                gene = gene_norm_resp.gene
                gene_norm_label = gene.label
                amplification_label = f"{gene_norm_label} Amplification"
                if all((sequence_id, start, end)):
                    # User provided input to make sequence location
                    seq_id, w = self.seqrepo_access.translate_identifier(
                        sequence_id, "ga4gh"
                    )
                    if w:
                        warnings.append(w)
                    else:
                        # Validate start/end are actually on the sequence
                        _, w = self.seqrepo_access.get_reference_sequence(
                            sequence_id, start=start, end=end
                        )
                        if w:
                            warnings.append(w)
                        else:
                            vrs_location = models.SequenceLocation(
                                sequenceReference=models.SequenceReference(
                                    refgetAccession=seq_id[0].split("ga4gh:")[-1]
                                ),
                                start=start - 1,
                                end=end,
                            )
                else:
                    # Use gene normalizer to get sequence location
                    seq_loc = get_priority_sequence_location(gene, self.seqrepo_access)
                    if seq_loc:
                        vrs_location = models.SequenceLocation(**seq_loc)
                    else:
                        warnings.append(
                            f"gene-normalizer could not find a priority sequence "
                            f"location for gene: {gene_norm_label}"
                        )

                if vrs_location:
                    vrs_location.id = ga4gh_identify(vrs_location)
                    vrs_cx = models.CopyNumberChange(
                        location=vrs_location,
                        copyChange=models.CopyChange.EFO_0030072.value,
                    )
                    vrs_cx.id = ga4gh_identify(vrs_cx)
                    variation = models.CopyNumberChange(
                        **vrs_cx.model_dump(exclude_none=True)
                    )
            else:
                warnings.append(f"gene-normalizer returned no match for gene: {gene}")

        return AmplificationToCxVarService(
            query=og_query,
            amplification_label=amplification_label,
            copy_number_change=variation,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
            ),
        )
