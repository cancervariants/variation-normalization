"""Module for general functionality throughout the app"""
from typing import List, Tuple, Optional, Dict
import re

from bioutils.sequences import aa1_to_aa3 as _aa1_to_aa3, aa3_to_aa1 as _aa3_to_aa1
from ga4gh.vrsatile.pydantic.vrs_models import Text
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor, GeneDescriptor
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models
from cool_seq_tool.data_sources import SeqRepoAccess

from variation.schemas.classification_response_schema import AmbiguousType
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.token_response_schema import Token, GnomadVcfToken
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.schemas.app_schemas import AmbiguousRegexType, Endpoint


def no_variation_entered() -> Tuple[None, List[str]]:
    """Return response when no variation queried.

    :return: None, list of warnings
    """
    warnings = ["No variation was entered to normalize"]
    return None, warnings


def no_variation_resp(
    label: str, _id: str, warnings: List, untranslatable_returns_text: bool = False
) -> Tuple[Optional[VariationDescriptor], List]:
    """Return Variation Descriptor with variation set as Text or return `None` for
    queries that could not be normalized

    :param str label: Initial input query
    :param str _id: _id field for variation descriptor
    :param List warnings: List of warnings
    :param bool untranslatable_returns_text: `True` return VRS Text Object when
        unable to translate or normalize query. `False` return `None` when
        unable to translate or normalize query.
    :return: Variation descriptor or `None`, warnings
    """
    if untranslatable_returns_text:
        text = models.Text(definition=label, type="Text")
        text._id = ga4gh_identify(text)
        variation = Text(**text.as_dict())
        resp = VariationDescriptor(id=_id, label=label, variation=variation)
    else:
        resp = None

    if not warnings:
        warnings.append(f"Unable to translate {label}")
    return resp, warnings


def _get_priority_sequence_location(locations: List[Dict],
                                    seqrepo_access: SeqRepoAccess) -> Optional[Dict]:
    """Get prioritized sequence location from list of locations
    Will prioritize GRCh8 over GRCh37. Will also only support chromosomes.

    :param List[Dict] locations: List of Chromosome and Sequence Locations represented
        as dictionaries
    :param SeqRepoAccess seqrepo_access: Client to access seqrepo
    :return: SequenceLocation represented as a dictionary if found
    """
    locs = [loc for loc in locations if loc["type"] == "SequenceLocation"]
    location = None
    if locs:
        if len(locs) > 1:
            loc38, loc37 = None, None
            for loc in locs:
                seq_id = loc["sequence_id"]
                aliases, _ = seqrepo_access.translate_identifier(seq_id)
                if aliases:
                    grch_aliases = [a for a in aliases
                                    if re.match(r"^GRCh3\d:chr(X|Y|\d+)$", a)]
                    if grch_aliases:
                        grch_alias = grch_aliases[0]
                        if grch_alias.startswith("GRCh38"):
                            loc38 = loc
                        elif grch_alias.startswith("GRCh37"):
                            loc37 = loc
                location = loc38 or loc37
        else:
            location = locs[0]

        if location:
            # DynamoDB stores as Decimal, so need to convert to int
            for k in {"start", "end"}:
                if location["interval"][k]["type"] == "Number":
                    location["interval"][k]["value"] = int(location["interval"][k]["value"])  # noqa: E501
    return location


def get_priority_sequence_location(gene_descriptor: GeneDescriptor,
                                   seqrepo_access: SeqRepoAccess) -> Optional[Dict]:
    """
    Get prioritized sequence location from a gene descriptor
    Will prioritize NCBI and then Ensembl. GRCh38 will be chosen over GRCh37.

    :param GeneDescriptor gene_descriptor: Gene descriptor
    :param SeqRepoAccess seqrepo_access: Client to access seqrepo
    :return: Prioritized sequence location represented as a dictionary if found
    """
    extensions = gene_descriptor.extensions or []

    # HGNC only provides ChromosomeLocation
    ensembl_loc, ncbi_loc = None, None
    for ext in extensions:
        if ext.name == "ncbi_locations":
            ncbi_loc = _get_priority_sequence_location(ext.value, seqrepo_access)
        elif ext.name == "ensembl_locations":
            ensembl_loc = _get_priority_sequence_location(ext.value, seqrepo_access)
    return ncbi_loc or ensembl_loc


def get_aa1_codes(aa: str) -> Optional[str]:
    """Get single AA codes given possible AA string (either single or three letter).
    Will also validate the input AA string

    :param aa: Input amino acid string
    :return: Amino acid string represented using single AA codes if valid. If invalid,
        will return None
    """
    aa1 = None
    if aa == "*":
        aa1 = aa
    else:
        try:
            # see if it's already 1 AA
            _aa1_to_aa3(aa)
        except KeyError:
            # see if it's 3 AA
            try:
                aa1 = _aa3_to_aa1(aa)
            except KeyError:
                pass
        else:
            aa1 = aa

    return aa1


def get_ambiguous_type(
    pos0, pos1, pos2, pos3, ambiguous_regex_type
) -> Optional[AmbiguousType]:
    ambiguous_type = None
    if ambiguous_regex_type == AmbiguousRegexType.REGEX_1:
        if all((
            isinstance(pos0, int),
            isinstance(pos1, int),
            isinstance(pos2, int),
            isinstance(pos3, int)
        )):
            ambiguous_type = AmbiguousType.AMBIGUOUS_1
        elif all((
            pos0 == "?",
            isinstance(pos1, int),
            isinstance(pos2, int),
            pos3 == "?"
        )):
            ambiguous_type = AmbiguousType.AMBIGUOUS_2
    elif ambiguous_regex_type == AmbiguousRegexType.REGEX_2:
        if all((
            pos0 == "?",
            isinstance(pos1, int),
            isinstance(pos2, int),
            pos3 is None
        )):
            ambiguous_type = AmbiguousType.AMBIGUOUS_5
    elif ambiguous_regex_type == AmbiguousRegexType.REGEX_3:
        if all((
            isinstance(pos0, int),
            pos1 is None,
            isinstance(pos2, int),
            pos3 == "?"
        )):
            ambiguous_type = AmbiguousType.AMBIGUOUS_7

    return ambiguous_type


def get_assembly(
    seqrepo_access, alt_ac: str
) -> Tuple[Optional[ClinVarAssembly], Optional[str]]:
    assembly = None
    warning = None

    grch38_aliases, _ = seqrepo_access.translate_identifier(alt_ac, "GRCh38")
    if grch38_aliases:
        assembly = ClinVarAssembly.GRCH38

    grch37_aliases, _ = seqrepo_access.translate_identifier(alt_ac, "GRCh37")
    if grch37_aliases:
        assembly = ClinVarAssembly.GRCH37

    if not assembly:
        warning = f"Unable to get GRCh37/GRCh38 assembly for: {alt_ac}"

    return assembly, warning


def get_hgvs_dup_del_mode(
    tokens: List[Token], endpoint_name: Endpoint,
    hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = None,
    baseline_copies: Optional[int] = None
) -> HGVSDupDelModeEnum:
    warning = None
    if len(tokens) == 1 and isinstance(tokens[0], GnomadVcfToken):
        # gnomad vcf should always be a literal seq expression (allele)
        hgvs_dup_del_mode = HGVSDupDelModeEnum.LITERAL_SEQ_EXPR
    else:
        if endpoint_name in {
            Endpoint.NORMALIZE, Endpoint.HGVS_TO_COPY_NUMBER_COUNT,
            Endpoint.HGVS_TO_COPY_NUMBER_CHANGE
        }:
            if hgvs_dup_del_mode:
                if hgvs_dup_del_mode == HGVSDupDelModeEnum.COPY_NUMBER_COUNT:
                    if not baseline_copies:
                        warning = f"{hgvs_dup_del_mode.value} mode requires `baseline_copies`"  # noqa: E501
                        return None, warning
            elif not hgvs_dup_del_mode and endpoint_name == Endpoint.NORMALIZE:
                hgvs_dup_del_mode = HGVSDupDelModeEnum.DEFAULT
            else:
                warning = (
                    f"hgvs_dup_del_mode must be either "
                    f"{HGVSDupDelModeEnum.COPY_NUMBER_COUNT.value} or "
                    f"{HGVSDupDelModeEnum.COPY_NUMBER_CHANGE.value}"
                )
                return None, warning
        else:
            hgvs_dup_del_mode = HGVSDupDelModeEnum.DEFAULT
    return hgvs_dup_del_mode, warning
