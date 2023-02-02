"""Module for to copy number variation translation"""
from typing import Tuple, Optional, List, Union
from datetime import datetime

from ga4gh.vrsatile.pydantic.vrs_models import AbsoluteCopyNumber, RelativeCopyNumber, \
    Text, RelativeCopyClass
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from pydantic import ValidationError
from gene.schemas import MatchType as GeneMatchType

from variation.to_vrs import ToVRS
from variation.schemas.app_schemas import Endpoint
from variation.schemas.hgvs_to_copy_number_schema import \
    HgvsToAbsoluteCopyNumberService, HgvsToRelativeCopyNumberService, \
    VALID_RELATIVE_COPY_CLASS
from variation.schemas.service_schema import AmplificationToRelCnvQuery, \
    AmplificationToRelCnvService, ClinVarAssembly, ParsedToAbsCnvQuery, \
    ParsedToAbsCnvService
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
        """Return response for invalid query for parsed_to_abs_cnv

        :param int start: Start position as residue coordinate
        :param int end: End position as residue coordinate
        :param int total_copies: Total copies for Absolute Copy Number variation object
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
    ) -> Tuple[Optional[Union[AbsoluteCopyNumber, RelativeCopyNumber, Text]], List[str]]:  # noqa: E501
        """Return copy number variation and warnings response

        :param HGVSDupDelModeEnum copy_number_type: The type of copy number variation.
            Must be either `absolute_cnv` or `relative_cnv`
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
            if copy_number_type == HGVSDupDelModeEnum.ABSOLUTE_CNV:
                variation = AbsoluteCopyNumber(**variation)
            else:
                variation = RelativeCopyNumber(**variation)
        return variation, warnings

    async def hgvs_to_absolute_copy_number(
        self, hgvs_expr: str, baseline_copies: int,
        do_liftover: bool = False, untranslatable_returns_text: bool = False
    ) -> HgvsToAbsoluteCopyNumberService:
        """Given hgvs, return abolute copy number variation

        :param str hgvs_expr: HGVS expression
        :param int baseline_copies: Baseline copies number
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: HgvsToAbsoluteCopyNumberService containing Absolute Copy Number
            Variation and warnings
        """
        validations, warnings = await self.get_validations(
            hgvs_expr, endpoint_name=Endpoint.HGVS_TO_ABSOLUTE_CN,
            hgvs_dup_del_mode=HGVSDupDelModeEnum.ABSOLUTE_CNV,
            baseline_copies=baseline_copies, do_liftover=do_liftover
        )
        abs_cnv, warnings = self._hgvs_to_cnv_resp(
            HGVSDupDelModeEnum.ABSOLUTE_CNV, hgvs_expr, do_liftover, validations,
            warnings, untranslatable_returns_text)

        return HgvsToAbsoluteCopyNumberService(
            hgvs_expr=hgvs_expr,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            ),
            absolute_copy_number=abs_cnv
        )

    async def hgvs_to_relative_copy_number(
        self, hgvs_expr: str, relative_copy_class: Optional[RelativeCopyClass],
        do_liftover: bool = False, untranslatable_returns_text: bool = False
    ) -> HgvsToRelativeCopyNumberService:
        """Given hgvs, return relative copy number variation

        :param str hgvs_expr: HGVS expression
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: HgvsToRelativeCopyNumberService containing Relative Copy Number
            Variation and warnings
        """
        if relative_copy_class and relative_copy_class.lower() not in VALID_RELATIVE_COPY_CLASS:  # noqa: E501
            return None, [f"{relative_copy_class} is not a valid relative copy class: "
                          f"{VALID_RELATIVE_COPY_CLASS}"]

        validations, warnings = await self.get_validations(
            hgvs_expr, endpoint_name=Endpoint.HGVS_TO_RELATIVE_CN,
            hgvs_dup_del_mode=HGVSDupDelModeEnum.RELATIVE_CNV,
            relative_copy_class=relative_copy_class, do_liftover=do_liftover
        )

        rel_cnv, warnings = self._hgvs_to_cnv_resp(
            HGVSDupDelModeEnum.RELATIVE_CNV, hgvs_expr, do_liftover, validations,
            warnings, untranslatable_returns_text)

        return HgvsToRelativeCopyNumberService(
            hgvs_expr=hgvs_expr,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            ),
            relative_copy_number=rel_cnv
        )

    def parsed_to_abs_cnv(
        self, start: int, end: int, total_copies: int,
        assembly: Optional[ClinVarAssembly] = None, chr: Optional[str] = None,
        accession: Optional[str] = None, untranslatable_returns_text: bool = False
    ) -> ParsedToAbsCnvService:
        """Given parsed ClinVar Copy Number Gain/Loss components, return Absolute
        Copy Number Variation

        :param int start: Start position as residue coordinate
        :param int end: End position as residue coordinate
        :param int total_copies: Total copies for Absolute Copy Number variation object
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
        :return: ParsedToAbsCnvService containing Absolute Copy Number variation
            and list of warnings
        """
        variation = None
        warnings = list()
        try:
            og_query = ParsedToAbsCnvQuery(
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
                    query = f"{assembly}:{chr}"
                    aliases, w = self.seqrepo_access.translate_identifier(query)
                    if w:
                        warnings.append(w)
                    else:
                        accession = ([a for a in aliases if a.startswith("refseq:")] or [None])[0]  # noqa: E501
                        if not accession:
                            warnings.append(f"Unable to find RefSeq accession for {query}")  # noqa: E501
                else:
                    warnings.append(f"{assembly} assembly is not current supported")
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
                    self.seqrepo_access.seqrepo_client[accession][start - 1]
                    self.seqrepo_access.seqrepo_client[accession][end]
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
                        "type": "AbsoluteCopyNumber",
                        "subject": location.as_dict(),
                        "copies": {"value": total_copies, "type": "Number"}
                    }
                    variation = self.hgvs_dup_del_mode._ga4gh_identify_cnv(
                        variation, is_abs=True)
                    variation = AbsoluteCopyNumber(**variation)

        return ParsedToAbsCnvService(
            query=og_query,
            absolute_copy_number=variation,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            )
        )

    def amplification_to_rel_cnv(
        self, gene: str, sequence_id: Optional[str] = None, start: Optional[int] = None,
        end: Optional[int] = None, untranslatable_returns_text: bool = False
    ) -> AmplificationToRelCnvService:
        """Return Relative Copy Number Variation for Amplification query
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
        :return: AmplificationToRelCnvService containing Relative Copy Number and
            list of warnings
        """
        warnings = list()
        amplification_label = None
        variation = None
        try:
            og_query = AmplificationToRelCnvQuery(
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
                        "type": "RelativeCopyNumber",
                        "subject": vrs_location.as_dict(),
                        "relative_copy_class": RelativeCopyClass.HIGH_LEVEL_GAIN.value
                    }
                    variation = self.hgvs_dup_del_mode._ga4gh_identify_cnv(
                        variation, is_abs=False)
                    variation = RelativeCopyNumber(**variation)
            else:
                warnings.append(f"gene-normalizer returned no match for gene: {gene}")

        if not variation and untranslatable_returns_text:
            text_variation = models.Text(definition=amplification_label, type="Text")
            text_variation.id = ga4gh_identify(text_variation)
            variation = Text(**text_variation.as_dict())

        return AmplificationToRelCnvService(
            query=og_query,
            amplification_label=amplification_label,
            relative_copy_number=variation,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()))
