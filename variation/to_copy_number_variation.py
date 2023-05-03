"""Module for to copy number variation translation"""
from typing import Tuple, Optional, List, Union
from datetime import datetime

from ga4gh.vrsatile.pydantic.vrs_models import CopyNumberCount, CopyNumberChange, \
    Text, CopyChange
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
    ParsedToCnVarService
from variation.schemas.validation_response_schema import ValidationSummary
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum, ServiceMeta
from variation.utils import get_mane_valid_result, get_priority_sequence_location
from variation.version import __version__


class ToCopyNumberVariation(ToVRS):
    """Class for representing copy number variation"""

    @staticmethod
    def _parsed_to_text(start: int, end: int, total_copies: int, warnings: List[str],
                        assembly: Optional[str] = None, chr: Optional[str] = None,
                        accession: Optional[str] = None) -> Tuple[Text, List[str]]:
        """Return response for invalid query for parsed_to_cn_var

        :param int start: Start position as residue coordinate
        :param int end: End position as residue coordinate
        :param int total_copies: Total copies for Copy Number Count variation object
        :param List[str] warnings: List of warnings
        :param Optional[ClinVarAssembly] assembly: Assembly. If `accession` is set,
            will ignore `assembly` and `chr`. If `accession` not set, must provide
            both `assembly` and `chr`.
        :param Optional[str] chr: Chromosome. Must set when `assembly` is set.
        :param Optional[str] accession: Accession. If `accession` is set,
            will ignore `assembly` and `chr`. If `accession` not set, must provide
            both `assembly` and `chr`.
        :return: Tuple containing text variation and warnings
        """
        text_label = ""
        for val, name in [(assembly, "assembly"), (chr, "chr"),
                          (accession, "accession"), (start, "start"),
                          (end, "end"), (total_copies, "total_copies")]:
            val = val if val else "None"
            text_label += f"{name}={val}&"

        definition = text_label[:-1]
        variation = models.Text(definition=definition, type="Text")
        _id = ga4gh_identify(variation)
        variation = Text(definition=definition, id=_id)
        return variation, warnings

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
                text._id = ga4gh_identify(text)
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

    def parsed_to_cn_var(
        self, start: int, end: int, total_copies: int,
        assembly: Optional[ClinVarAssembly] = None, chr: Optional[str] = None,
        accession: Optional[str] = None, untranslatable_returns_text: bool = False
    ) -> ParsedToCnVarService:
        """Given parsed ClinVar Copy Number Gain/Loss components, return Copy Number
        Count Variation

        :param int start: Start position as residue coordinate
        :param int end: End position as residue coordinate
        :param int total_copies: Total copies for Copy Number Count variation object
        :param Optional[ClinVarAssembly] assembly: Assembly. If `accession` is set,
            will ignore `assembly` and `chr`. If `accession` not set, must provide
            both `assembly` and `chr`.
        :param Optional[str] chr: Chromosome. Must set when `assembly` is set.
        :param Optional[str] accession: Accession. If `accession` is set,
            will ignore `assembly` and `chr`. If `accession` not set, must provide
            both `assembly` and `chr`.
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: ParsedToCnVarService containing Copy Number Count variation
            and list of warnings
        """
        variation = None
        warnings = list()
        try:
            og_query = ParsedToCnVarQuery(
                assembly=assembly, chr=chr, accession=accession, start=start, end=end,
                total_copies=total_copies)
        except ValidationError as e:
            warnings.append(str(e))
            og_query = None
        else:
            if accession:
                pass
            elif assembly and chr:
                if assembly == ClinVarAssembly.HG38:
                    assembly = ClinVarAssembly.GRCH38
                elif assembly == ClinVarAssembly.HG19:
                    assembly = ClinVarAssembly.GRCH37
                elif assembly == ClinVarAssembly.HG18:
                    assembly = ClinVarAssembly.NCBI36

                if assembly != ClinVarAssembly.NCBI36:
                    # Variation Normalizer does not support NCBI36 yet
                    query = f"{assembly.value}:{chr}"
                    aliases, w = self.seqrepo_access.translate_identifier(query)
                    if w:
                        warnings.append(w)
                    else:
                        accession = ([a for a in aliases if a.startswith("refseq:")] or [None])[0]  # noqa: E501
                        if not accession:
                            warnings.append(f"Unable to find RefSeq accession for {query}")  # noqa: E501
                else:
                    warnings.append(
                        f"{assembly.value} assembly is not currently supported"
                    )
            else:
                warnings.append("Must provide either `accession` or both `assembly` "
                                "and `chr`.")

        if warnings:
            if untranslatable_returns_text:
                variation, warnings = self._parsed_to_text(
                    start, end, total_copies, warnings, assembly, chr, accession)
        else:
            try:
                sequence_id, w = self.seqrepo_access.translate_identifier(accession,
                                                                          "ga4gh")
            except (IndexError, TypeError):
                warnings.append(f"{accession} does not have an associated "
                                f"ga4gh identifier")
            else:
                if w:
                    warnings.append(w)
                else:
                    sequence_id = sequence_id[0]

            if warnings:
                if untranslatable_returns_text:
                    variation, warnings = self._parsed_to_text(
                        start, end, total_copies, warnings, assembly, chr, accession)
            else:
                try:
                    self.seqrepo_access.sr[accession][start - 1]
                    self.seqrepo_access.sr[accession][end]
                except ValueError as e:
                    warnings.append(str(e).replace("start", "Position"))
                    if untranslatable_returns_text:
                        variation, warnings = self._parsed_to_text(
                            start, end, total_copies, warnings, assembly, chr,
                            accession)
                else:
                    location = models.SequenceLocation(
                        type="SequenceLocation",
                        sequence_id=sequence_id,
                        interval=models.SequenceInterval(
                            type="SequenceInterval",
                            start=models.IndefiniteRange(
                                comparator="<=",
                                value=start - 1,
                                type="IndefiniteRange"),
                            end=models.IndefiniteRange(
                                comparator=">=",
                                value=end,
                                type="IndefiniteRange")
                        )
                    )
                    location._id = ga4gh_identify(location)
                    variation = {
                        "type": "CopyNumberCount",
                        "subject": location.as_dict(),
                        "copies": {"value": total_copies, "type": "Number"}
                    }
                    variation["_id"] = ga4gh_identify(
                        models.CopyNumberCount(**variation)
                    )
                    variation = CopyNumberCount(**variation)

        return ParsedToCnVarService(
            query=og_query,
            copy_number_count=variation,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            )
        )

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
                                type="SequenceLocation",
                                sequence_id=seq_id[0],
                                interval=models.SequenceInterval(
                                    type="SequenceInterval",
                                    start=models.Number(type="Number", value=start - 1),
                                    end=models.Number(type="Number", value=end)))
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
                    vrs_location._id = ga4gh_identify(vrs_location)
                    variation = {
                        "type": "CopyNumberChange",
                        "subject": vrs_location.as_dict(),
                        "copy_change": CopyChange.HIGH_LEVEL_GAIN.value
                    }
                    variation["_id"] = ga4gh_identify(
                        models.CopyNumberChange(**variation)
                    )
                    variation = CopyNumberChange(**variation)
            else:
                warnings.append(f"gene-normalizer returned no match for gene: {gene}")

        if not variation and untranslatable_returns_text:
            text_variation = models.Text(definition=amplification_label, type="Text")
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
