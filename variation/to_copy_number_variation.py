"""Module for to copy number variation translation"""
from typing import Tuple, Optional, List, Union, Dict, Literal
from datetime import datetime

from ga4gh.vrsatile.pydantic.vrs_models import CopyNumberCount, CopyNumberChange, \
    Text, CopyChange, VRSTypes
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from pydantic import ValidationError
from gene.schemas import MatchType as GeneMatchType

from variation.to_vrs import ToVRS
from variation.schemas.app_schemas import Endpoint
from variation.schemas.hgvs_to_copy_number_schema import \
    HgvsToCopyNumberCountService, HgvsToCopyNumberChangeService, \
    VALID_COPY_CHANGE
from variation.schemas.service_schema import AmplificationToCxVarQuery, \
    AmplificationToCxVarService, ClinVarAssembly, ParsedToCnVarQuery, \
    ParsedToCnVarService, ParsedToCxVarQuery, ParsedToCxVarService
from variation.schemas.validation_response_schema import ValidationSummary
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum, ServiceMeta
from variation.utils import get_mane_valid_result, get_priority_sequence_location
from variation.version import __version__


class ToCopyNumberVariation(ToVRS):
    """Class for representing copy number variation"""

    @staticmethod
    def _parsed_to_text(
        params: Dict
    ) -> Text:
        """Return response for invalid query for parsed_to_cn_var

        :param params: Parameters for initial query
        :return: Variation represented as VRS Text object
        """
        text_label = ""
        for name, val in params.items():
            val = val if val else "None"
            text_label += f"{name}={val}&"

        definition = text_label[:-1]
        variation = models.Text(definition=definition, type="Text")
        _id = ga4gh_identify(variation)
        variation = Text(definition=definition, id=_id)
        return variation

    def _hgvs_to_cnv_resp(
        self, copy_number_type: HGVSDupDelModeEnum, hgvs_expr: str, do_liftover: bool,
        validations: Tuple[Optional[ValidationSummary], Optional[List[str]]],
        warnings: List[str], untranslatable_returns_text: bool = False
    ) -> Tuple[Optional[Union[CopyNumberCount, CopyNumberChange, Text]], List[str]]:  # noqa: E501
        """Return copy number variation and warnings response

        :param HGVSDupDelModeEnum copy_number_type: The type of copy number variation.
            Must be either `copy_number_count` or `copy_number_change`
        :param str hgvs_expr: HGVS expression
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :param Tuple[Optional[ValidationSummary], Optional[List[str]]]: Validation
            summary and warnings for hgvs_expr
        :param List[str] warnings: List of warnings
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: CopyNumberVariation and warnings
        """
        variation = None
        if do_liftover:
            valid_result = get_mane_valid_result(hgvs_expr, validations, [])
            if valid_result:
                variation = valid_result.variation
            else:
                warnings.append(f"Unable to translate {hgvs_expr} to "
                                f"copy number variation")
        else:
            translations, warnings = self.get_translations(validations, warnings)
            if translations:
                variation = translations[0]

        if not variation:
            if hgvs_expr and hgvs_expr.strip() and untranslatable_returns_text:
                text = models.Text(definition=hgvs_expr, type="Text")
                text.id = ga4gh_identify(text)
                variation = Text(**text.as_dict())
        else:
            if copy_number_type == HGVSDupDelModeEnum.COPY_NUMBER_COUNT:
                variation = CopyNumberCount(**variation)
            else:
                variation = CopyNumberChange(**variation)
        return variation, warnings

    async def hgvs_to_copy_number_count(
        self, hgvs_expr: str, baseline_copies: int,
        do_liftover: bool = False, untranslatable_returns_text: bool = False
    ) -> HgvsToCopyNumberCountService:
        """Given hgvs, return abolute copy number variation

        :param str hgvs_expr: HGVS expression
        :param int baseline_copies: Baseline copies number
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: HgvsToCopyNumberCountService containing Copy Number Count
            Variation and warnings
        """
        validations, warnings = await self.get_validations(
            hgvs_expr, endpoint_name=Endpoint.HGVS_TO_COPY_NUMBER_COUNT,
            hgvs_dup_del_mode=HGVSDupDelModeEnum.COPY_NUMBER_COUNT,
            baseline_copies=baseline_copies, do_liftover=do_liftover
        )
        cn_var, warnings = self._hgvs_to_cnv_resp(
            HGVSDupDelModeEnum.COPY_NUMBER_COUNT, hgvs_expr, do_liftover, validations,
            warnings, untranslatable_returns_text)

        return HgvsToCopyNumberCountService(
            hgvs_expr=hgvs_expr,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            ),
            copy_number_count=cn_var
        )

    async def hgvs_to_copy_number_change(
        self, hgvs_expr: str, copy_change: Optional[CopyChange],
        do_liftover: bool = False, untranslatable_returns_text: bool = False
    ) -> HgvsToCopyNumberChangeService:
        """Given hgvs, return copy number change variation

        :param str hgvs_expr: HGVS expression
        :param Optional[CopyChange] copy_change: The copy change
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: HgvsToCopyNumberChangeService containing Copy Number Change
            Variation and warnings
        """
        if copy_change and copy_change.lower() not in VALID_COPY_CHANGE:  # noqa: E501
            return None, [f"{copy_change} is not a valid copy change: "
                          f"{VALID_COPY_CHANGE}"]

        validations, warnings = await self.get_validations(
            hgvs_expr, endpoint_name=Endpoint.HGVS_TO_COPY_NUMBER_CHANGE,
            hgvs_dup_del_mode=HGVSDupDelModeEnum.COPY_NUMBER_CHANGE,
            copy_change=copy_change, do_liftover=do_liftover
        )

        cx_var, warnings = self._hgvs_to_cnv_resp(
            HGVSDupDelModeEnum.COPY_NUMBER_CHANGE, hgvs_expr, do_liftover, validations,
            warnings, untranslatable_returns_text)

        return HgvsToCopyNumberChangeService(
            hgvs_expr=hgvs_expr,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            ),
            copy_number_change=cx_var
        )

    def _get_parsed_ac(
        self, assembly: ClinVarAssembly, chr: str
    ) -> Tuple[Optional[str], Optional[str]]:
        """Get accession for parsed components

        :param assembly: Assembly
        :param chr: Chromosome
        :return: Tuple containing accession (if successful) and warnings (if error)
        """
        accession = None
        warning = None

        if assembly == ClinVarAssembly.HG38:
            assembly = ClinVarAssembly.GRCH38
        elif assembly == ClinVarAssembly.HG19:
            assembly = ClinVarAssembly.GRCH37
        elif assembly == ClinVarAssembly.HG18:
            assembly = ClinVarAssembly.NCBI36

        if assembly != ClinVarAssembly.NCBI36:
            # Variation Normalizer does not support NCBI36 yet
            query = f"{assembly.value}:{chr}"
            aliases, warning = self.seqrepo_access.translate_identifier(query)
            if not warning:
                accession = ([a for a in aliases if a.startswith("refseq:")] or [None])[0]  # noqa: E501
                if not accession:
                    warning = f"Unable to find RefSeq accession for {query}"
        else:
            warning = f"{assembly.value} assembly is not currently supported"

        return accession, warning

    def _validate_pos(self, accession: str, pos: int) -> Optional[str]:
        """Validate position for parsed components

        :param accession: Genomic accession
        :param pos: Position on accession
        :return: Warning if invalid position or sequence
        """
        try:
            self.seqrepo_access.sr[accession][pos - 1]
        except ValueError as e:
            warning = f"SeqRepo ValueError: {str(e).replace('start', 'Position')}"
        else:
            warning = None
        return warning

    def _get_vrs_loc_start_or_end(
        self, accession: str, pos0: int,
        pos_type: Literal[
            VRSTypes.NUMBER, VRSTypes.DEFINITE_RANGE, VRSTypes.INDEFINITE_RANGE
        ],
        is_start: bool = True, pos1: Optional[int] = None
    ) -> Tuple[
        Optional[Union[models.Number, models.DefiniteRange, models.IndefiniteRange]],
        Optional[str]
    ]:
        """Get VRS Sequence Location start and end values

        :param accession: Genomic accession for sequence
        :param pos0: Position (residue coords). If `pos_type` is a definite range,
            this will be the min start position
        :param pos_type: Type of the pos value in VRS Sequence Location
        :param is_start: `True` if position(s) describing VRS start value. `False` if
            position(s) describing VRS end value
        :param pos1: Only set when end is a definite range, this will be the max end
            position
        :return: Tuple containing VRS start or end value for sequence location
            (if valid) and warning (if invalid)
        """
        vrs_val = None
        warning = None

        if pos_type == VRSTypes.DEFINITE_RANGE:
            warning = self._validate_pos(accession, pos1)
            if warning:
                return vrs_val, warning

            vrs_val = models.DefiniteRange(
                min=pos0 - 1 if is_start else pos0 + 1,
                max=pos1 - 1 if is_start else pos1 + 1,
                type="DefiniteRange"
            )
        elif pos_type == VRSTypes.NUMBER:
            vrs_val = models.Number(
                value=pos0 - 1 if is_start else pos0,
                type="Number"
            )
        else:
            vrs_val = models.IndefiniteRange(
                comparator="<=" if is_start else ">=",
                value=pos0 - 1 if is_start else pos0,
                type="IndefiniteRange"
            )

        return vrs_val, warning

    def _get_parsed_seq_loc(
        self, accession: str, start0: int,
        start_pos_type: Union[
            VRSTypes.NUMBER, VRSTypes.DEFINITE_RANGE, VRSTypes.INDEFINITE_RANGE
        ],
        end0: int,
        end_pos_type: Union[
            VRSTypes.NUMBER, VRSTypes.DEFINITE_RANGE, VRSTypes.INDEFINITE_RANGE
        ],
        start1: Optional[int] = None, end1: Optional[int] = None
    ) -> Tuple[Optional[Dict], Optional[str]]:
        """Get sequence location for parsed components

        :param accession: Genomic accession for sequence
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
        :return: Tuple containing VRS sequence location represented as dict (if valid)
            and warning (if invalid)
        """
        seq_loc = None
        warning = None

        try:
            sequence_id, warning = self.seqrepo_access.translate_identifier(
                accession, "ga4gh"
            )
        except (IndexError, TypeError):
            warning = f"{accession} does not have an associated ga4gh identifier"
        else:
            if not warning:
                sequence_id = sequence_id[0]

                for pos in [start0, end0]:
                    # validate start0 and end0 since they're always required
                    warning = self._validate_pos(accession, pos)
                    if warning:
                        return seq_loc, warning

                start_vrs, warning = self._get_vrs_loc_start_or_end(
                    accession, start0, start_pos_type, is_start=True, pos1=start1
                )
                if warning:
                    return seq_loc, warning

                end_vrs, warning = self._get_vrs_loc_start_or_end(
                    accession, end0, end_pos_type, is_start=False, pos1=end1
                )
                if warning:
                    return seq_loc, warning

                seq_loc = models.SequenceLocation(
                    type="SequenceLocation",
                    sequence_id=sequence_id,
                    start=start_vrs,
                    end=end_vrs
                )
                seq_loc.id = ga4gh_identify(seq_loc)

        return seq_loc.as_dict() if seq_loc else seq_loc, warning

    def parsed_to_copy_number(
        self,
        start0: int,
        end0: int,
        copy_number_type: Literal[
            VRSTypes.COPY_NUMBER_CHANGE, VRSTypes.COPY_NUMBER_COUNT
        ],
        total_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
        start_pos_type: Literal[
            VRSTypes.NUMBER, VRSTypes.DEFINITE_RANGE, VRSTypes.INDEFINITE_RANGE
        ] = VRSTypes.NUMBER,
        end_pos_type: Literal[
            VRSTypes.NUMBER, VRSTypes.DEFINITE_RANGE, VRSTypes.INDEFINITE_RANGE
        ] = VRSTypes.NUMBER,
        start1: Optional[int] = None,
        end1: Optional[int] = None,
        assembly: Optional[ClinVarAssembly] = None,
        chr: Optional[str] = None,
        accession: Optional[str] = None,
        untranslatable_returns_text: bool = False
    ) -> Union[ParsedToCnVarService, ParsedToCxVarService]:
        """Given parsed genomic components, return Copy Number Count or Copy Number
        Change Variation

        :param start0: Start position (residue coords). If start is a definite range,
            this will be the min start position
        :param end0: End position (residue coords). If end is a definite range, this
            will be the min end position
        :param copy_number_type: The kind of Copy Number you wish to represent
        :param total_copies: Total copies for Copy Number Count variation object.
            Only used when `copy_number_type` is Copy Number Count.
        :param copy_change: Copy Change.
            Only used when `copy_number_type` is Copy Number Change.
        :param assembly: Assembly. Ignored, along with `chr`, if `accession` is set.
            If `accession` not set, must provide both `assembly` and `chr`.
        :param chr: Chromosome. Must set when `assembly` is set.
        :param accession: Genomic accession. If `accession` is set, will ignore
            `assembly` and `chr`. If `accession` not set, must provide both `assembly`
            and `chr`.
        :param start_pos_type: Type of the start value in VRS Sequence Location
        :param end_pos_type: Type of the end value in VRS Sequence Location
        :param start1: Only set when start is a definite range, this will be the max
            start position
        :param end1: Only set when end is a definite range, this will be the max end
            position
        :param untranslatable_returns_text: `True` return VRS Text Object when unable to
            translate or normalize query. `False` return `None` when unable to translate
            or normalize query.
        :return: If `copy_number_type` is Copy Number Count, return ParsedToCnVarService
            containing Copy Number Count variation and list of warnings. Else, return
            ParsedToCxVarService containing Copy Number Change variation and list of
            warnings
        """
        variation = None
        warnings = []
        params = {
            "assembly": assembly,
            "chr": chr,
            "accession": accession,
            "start0": start0,
            "end0": end0,
            "start_pos_type": start_pos_type,
            "end_pos_type": end_pos_type,
            "start1": start1,
            "end1": end1
        }

        is_cx = copy_number_type == VRSTypes.COPY_NUMBER_CHANGE
        if is_cx:
            params["copy_change"] = copy_change
        else:
            params["total_copies"] = total_copies

        try:
            og_query = ParsedToCxVarQuery(**params) if is_cx else ParsedToCnVarQuery(**params)  # noqa: E501
        except ValidationError as e:
            warnings.append(str(e))
            og_query = None
        else:
            if not accession:
                accession, warning = self._get_parsed_ac(assembly, chr)
            else:
                warning = None

            if warning:
                warnings.append(warning)
            else:
                seq_loc, warning = self._get_parsed_seq_loc(
                    accession, start0, start_pos_type, end0, end_pos_type,
                    start1=start1, end1=end1
                )
                if warning:
                    warnings.append(warning)

        if warnings:
            if untranslatable_returns_text:
                variation = self._parsed_to_text(params)
        else:
            if is_cx:
                variation = {
                    "type": "CopyNumberChange",
                    "subject": seq_loc,
                    "copy_change": copy_change
                }
                variation["id"] = ga4gh_identify(models.CopyNumberChange(**variation))
                variation = CopyNumberChange(**variation)
            else:
                variation = {
                    "type": "CopyNumberCount",
                    "subject": seq_loc,
                    "copies": {"value": total_copies, "type": "Number"}
                }
                variation["id"] = ga4gh_identify(models.CopyNumberCount(**variation))
                variation = CopyNumberCount(**variation)

        service_params = {
            "og_query": og_query,
            "warnings": warnings,
            "service_meta_": ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            )
        }

        if is_cx:
            service_params["copy_number_change"] = variation
        else:
            service_params["copy_number_count"] = variation

        return ParsedToCxVarService(**service_params) if is_cx else ParsedToCnVarService(**service_params)  # noqa: E501

    def amplification_to_cx_var(
        self, gene: str, sequence_id: Optional[str] = None, start: Optional[int] = None,
        end: Optional[int] = None, untranslatable_returns_text: bool = False
    ) -> AmplificationToCxVarService:
        """Return Copy Number Change Variation for Amplification query
        Parameter priority:
            1. sequence_id, start, end (must provide ALL)
            2. use the gene-normalizer to get the SequenceLocation

        :param str gene: Gene query
        :param Optional[str] sequence_id: Sequence ID for the location. If set,
            must also provide `start` and `end`
        :param Optional[int] start: Start position as residue coordinate for the
            sequence location. If set, must also provide `sequence_id` and `end`
        :param Optional[int] end: End position as residue coordinate for the sequence
            location. If set, must also provide `sequence_id` and `start`
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: AmplificationToCxVarService containing Copy Number Change and
            list of warnings
        """
        warnings = list()
        amplification_label = None
        variation = None
        try:
            og_query = AmplificationToCxVarQuery(
                gene=gene, sequence_id=sequence_id, start=start, end=end)
        except ValidationError as e:
            warnings.append(str(e))
            og_query = None
        else:
            # Need to validate the input gene
            gene_norm_resp = self.gene_normalizer.normalize(gene)
            if gene_norm_resp.match_type != GeneMatchType.NO_MATCH:
                vrs_location = None
                gene_descriptor = gene_norm_resp.gene_descriptor
                gene_norm_label = gene_descriptor.label
                amplification_label = f"{gene_norm_label} Amplification"
                if all((sequence_id, start, end)):
                    # User provided input to make sequence location
                    seq_id, w = self.seqrepo_access.translate_identifier(
                        sequence_id, "ga4gh")
                    if w:
                        warnings.append(w)
                    else:
                        # Validate start/end are actually on the sequence
                        _, w = self.seqrepo_access.get_reference_sequence(
                            sequence_id, start, end)
                        if w:
                            warnings.append(w)
                        else:
                            vrs_location = models.SequenceLocation(
                                sequence_id=seq_id[0],
                                start=models.Number(value=start - 1),
                                end=models.Number(value=end))
                else:
                    # Use gene normalizer to get sequence location
                    seq_loc = get_priority_sequence_location(
                        gene_descriptor, self.seqrepo_access)
                    if seq_loc:
                        vrs_location = models.SequenceLocation(**seq_loc)
                    else:
                        warnings.append(
                            f"gene-normalizer could not find a priority sequence "
                            f"location for gene: {gene_norm_label}")

                if vrs_location:
                    vrs_location.id = ga4gh_identify(vrs_location)
                    vrs_cx = models.CopyNumberChange(
                        subject=vrs_location,
                        copy_change=CopyChange.HIGH_LEVEL_GAIN.value)
                    vrs_cx.id = ga4gh_identify(vrs_cx)
                    variation = CopyNumberChange(**vrs_cx.as_dict())
            else:
                warnings.append(f"gene-normalizer returned no match for gene: {gene}")

        if not variation and untranslatable_returns_text:
            text_variation = models.Text(definition=amplification_label)
            text_variation.id = ga4gh_identify(text_variation)
            variation = Text(**text_variation.as_dict())

        return AmplificationToCxVarService(
            query=og_query,
            amplification_label=amplification_label,
            copy_number_change=variation,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()))
