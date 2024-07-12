"""Module for general functionality throughout the app"""

import contextlib
import re
from typing import Dict, List, Literal, Optional, Tuple, Union

from bioutils.sequences import aa1_to_aa3 as _aa1_to_aa3
from bioutils.sequences import aa3_to_aa1 as _aa3_to_aa1
from cool_seq_tool.handlers import SeqRepoAccess
from ga4gh.core import core_models

from variation.schemas.app_schemas import AmbiguousRegexType
from variation.schemas.classification_response_schema import AmbiguousType
from variation.schemas.service_schema import ClinVarAssembly


def update_warnings_for_no_resp(label: str, warnings: List[str]) -> None:
    """Mutate `warnings` when unable to return a response

    :param label: Initial input query
    :param warnings: List of warnings to mutate
    """
    if not warnings:
        warnings.append(f"Unable to translate {label}")


def _get_priority_sequence_location(
    locations: List[Dict], seqrepo_access: SeqRepoAccess
) -> Optional[Dict]:
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
                seq_id = f"ga4gh:{loc['sequenceReference']['refgetAccession']}"
                aliases, _ = seqrepo_access.translate_identifier(seq_id)
                if aliases:
                    grch_aliases = [
                        a for a in aliases if re.match(r"^GRCh3\d:chr(X|Y|\d+)$", a)
                    ]
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
                location[k] = int(location[k])
    return location


def get_priority_sequence_location(
    gene: core_models.Gene, seqrepo_access: SeqRepoAccess
) -> Optional[Dict]:
    """Get prioritized sequence location from a gene
    Will prioritize NCBI and then Ensembl. GRCh38 will be chosen over GRCh37.

    :param gene: GA4GH Core Gene
    :param seqrepo_access: Client to access seqrepo
    :return: Prioritized sequence location represented as a dictionary if found
    """
    extensions = gene.extensions or []

    # HGNC only provides ChromosomeLocation
    ensembl_loc, ncbi_loc = None, None
    for ext in extensions:
        if ext.name == "ncbi_locations":
            ncbi_loc = _get_priority_sequence_location(ext.value, seqrepo_access)
        elif ext.name == "ensembl_locations":
            ensembl_loc = _get_priority_sequence_location(ext.value, seqrepo_access)
    return ncbi_loc or ensembl_loc


def get_aa1_codes(aa: str) -> Optional[str]:
    """Get 1 letter AA codes given possible AA string (either 1 or 3 letter).
    Will also validate the input AA string.

    :param aa: Input amino acid string. Case sensitive.
    :return: Amino acid string represented using 1 letter AA codes if valid. If invalid,
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
            with contextlib.suppress(KeyError):
                aa1 = _aa3_to_aa1(aa)
        else:
            aa1 = aa

    return aa1


def get_ambiguous_type(
    pos0: Union[int, Literal["?"]],
    pos1: Optional[Union[int, Literal["?"]]],
    pos2: Union[int, Literal["?"]],
    pos3: Optional[Union[int, Literal["?"]]],
    ambiguous_regex_type: AmbiguousRegexType,
) -> Optional[AmbiguousType]:
    """Get the ambiguous type given positions and regex used

    :param pos0: Position 0
    :param pos1: Position 1
    :param pos2: Position 2
    :param pos3: Position 3
    :param ambiguous_regex_type: The matched regex pattern
    :return: Corresponding ambiguous type if a match is found. Else, `None`
    """
    ambiguous_type = None
    if ambiguous_regex_type == AmbiguousRegexType.REGEX_1:
        if all(
            (
                isinstance(pos0, int),
                isinstance(pos1, int),
                isinstance(pos2, int),
                isinstance(pos3, int),
            )
        ):
            ambiguous_type = AmbiguousType.AMBIGUOUS_1
        elif all(
            (pos0 == "?", isinstance(pos1, int), isinstance(pos2, int), pos3 == "?")
        ):
            ambiguous_type = AmbiguousType.AMBIGUOUS_2
    elif ambiguous_regex_type == AmbiguousRegexType.REGEX_2:
        if all(
            (pos0 == "?", isinstance(pos1, int), isinstance(pos2, int), pos3 is None)
        ):
            ambiguous_type = AmbiguousType.AMBIGUOUS_5
    elif all(
        (
            ambiguous_regex_type == AmbiguousRegexType.REGEX_3,
            isinstance(pos0, int),
            pos1 is None,
            isinstance(pos2, int),
            pos3 == "?",
        )
    ):
        ambiguous_type = AmbiguousType.AMBIGUOUS_7

    return ambiguous_type


def get_assembly(
    seqrepo_access: SeqRepoAccess, alt_ac: str
) -> Tuple[Optional[ClinVarAssembly], Optional[str]]:
    """Get GRCh assembly for given genomic RefSeq accession

    :param seqrepo_access: Access to SeqRepo client
    :param alt_ac: Genomic RefSeq accession
    :return: Tuple containing the corresponding GRCh assembly, if found in SeqRepo and
        optional warning message if an assembly is not found in SeqRepo
    """
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


def get_refget_accession(
    seqrepo_access: SeqRepoAccess, alias: str, errors: List[str]
) -> Optional[str]:
    """Get refget accession for a given alias

    :param seqrepo_access: Access to SeqRepo client
    :param alias: Alias to translate
    :param errors: List of errors. This will get mutated if an error occurs.
    :return: RefGet Accession if successful, else `None`
    """
    refget_accession = None
    try:
        ids = seqrepo_access.translate_sequence_identifier(alias, "ga4gh")
    except KeyError as e:
        errors.append(str(e))
    else:
        if not ids:
            errors.append(f"Unable to find ga4gh sequence identifiers for: {alias}")

        refget_accession = ids[0].split("ga4gh:")[-1]
    return refget_accession
