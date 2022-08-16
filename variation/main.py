"""Main application for FastAPI."""
from enum import Enum
from typing import Dict, List, Optional, Union
from datetime import datetime
from urllib.parse import unquote

import pkg_resources
from fastapi import FastAPI, Query
from fastapi.openapi.utils import get_openapi
from pydantic import ValidationError
import python_jsonschema_objects
from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass
from hgvs.exceptions import HGVSError
from bioutils.exceptions import BioutilsError
from ga4gh.vrs import models

from variation.schemas import ToVRSService, NormalizeService, ServiceMeta
from variation.schemas.hgvs_to_copy_number_schema import \
    HgvsToAbsoluteCopyNumberService, HgvsToRelativeCopyNumberService
from variation.query import QueryHandler
from variation.schemas.normalize_response_schema \
    import HGVSDupDelMode as HGVSDupDelModeEnum, ToCanonicalVariationFmt, \
    ToCanonicalVariationService, TranslateIdentifierService
from variation.schemas.service_schema import ClinVarAssembly, ParsedToAbsCnvService
from .version import __version__
from .schemas.vrs_python_translator_schema import TranslateFromFormat, \
    TranslateFromService, TranslateFromQuery, TranslateToHGVSQuery, TranslateToQuery,\
    TranslateToService, VrsPythonMeta


class Tags(Enum):
    """Define tags for endpoints"""

    SEQREPO = "SeqRepo"
    TO_PROTEIN_VARIATION = "To Protein Variation"
    TO_CANONICAL = "To Canonical Variation"
    VRS_PYTHON = "VRS-Python"
    TO_COPY_NUMBER_VARIATION = "To Copy Number Variation"


app = FastAPI(
    docs_url="/variation",
    openapi_url="/variation/openapi.json",
    swagger_ui_parameters={"tryItOutEnabled": True}
)
query_handler = QueryHandler()


def custom_openapi() -> Dict:
    """Generate custom fields for OpenAPI response."""
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="The VICC Variation Normalizer",
        version=__version__,
        description="Services and guidelines for normalizing variations.",
        routes=app.routes
    )

    openapi_schema["info"]["contact"] = {
        "name": "Alex H. Wagner",
        "email": "Alex.Wagner@nationwidechildrens.org",
        "url": "https://www.nationwidechildrens.org/specialties/institute-for-genomic-medicine/research-labs/wagner-lab"  # noqa: E501
    }
    app.openapi_schema = openapi_schema
    return app.openapi_schema


app.openapi = custom_openapi

to_vrs_untranslatable_descr = ("`True` returns VRS Text object when unable to translate"
                               " or normalize query. `False` returns an empty list "
                               "when unable to translate or normalize query.")
translate_summary = ("Translate a human readable variation description to VRS"
                     " variation(s).")
translate_description = ("Translate a human readable variation description to "
                         "VRS variation(s)."
                         " Performs fully-justified allele normalization. "
                         " Does not do any liftover operations or make any inferences "
                         "about the query.")
translate_response_description = "A  response to a validly-formed query."
q_description = "Human readable variation description on GRCh37 or GRCh38 assembly"


@app.get("/variation/to_vrs",
         summary=translate_summary,
         response_description=translate_response_description,
         response_model=ToVRSService,
         description=translate_description
         )
async def to_vrs(
    q: str = Query(..., description=q_description),
    untranslatable_returns_text: bool = Query(False,
                                              description=to_vrs_untranslatable_descr)
) -> ToVRSService:
    """Translate a human readable variation description to VRS variation(s).
        Performs fully-justified allele normalization. Does not do any liftover
        operations or make any inferences about the query.

    :param str q: Human readable variation description on GRCh37 or GRCh38 assembly
    :param bool untranslatable_returns_text: `True` return VRS Text Object when
        unable to translate or normalize query. `False` returns empty list when
        unable to translate or normalize query.
    :return: ToVRSService model for variation
    """
    resp = await query_handler.to_vrs_handler.to_vrs(unquote(q),
                                                     untranslatable_returns_text)
    return resp

untranslatable_descr = ("`True` returns VRS Text object when unable to translate or "
                        "normalize query. `False` returns `None` when unable to "
                        "translate or normalize query.")

normalize_summary = ("Normalizes and translates a human readable variation description "
                     "to a single VRSATILE Variation Descriptor.")
normalize_response_description = "A response to a validly-formed query."
normalize_description = ("Normalizes and translates a human readable variation "
                         "description to a single VRSATILE Variation Descriptor. "
                         "Performs fully-justified allele normalization. Will liftover"
                         " to GRCh38 and aligns to a priority transcript. Will make "
                         "inferences about the query.")
q_description = "Human readable variation description on GRCh37 or GRCh38 assembly"
hgvs_dup_del_mode_decsr = ("Must be one of: `default`, `absolute_cnv`, `relative_cnv`, "
                           "`repeated_seq_expr`, `literal_seq_expr`. This"
                           " parameter determines how to interpret HGVS "
                           "dup/del expressions in VRS.")


@app.get("/variation/normalize",
         summary=normalize_summary,
         response_description=normalize_response_description,
         response_model=NormalizeService,
         description=normalize_description
         )
async def normalize(
    q: str = Query(..., description=q_description),
    hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = Query(
        HGVSDupDelModeEnum.DEFAULT, description=hgvs_dup_del_mode_decsr),
    baseline_copies: Optional[int] = Query(
        None, description="Baseline copies for HGVS duplications and deletions represented as Absolute Copy Number Variation"),  # noqa: E501
    relative_copy_class: Optional[RelativeCopyClass] = Query(
        None, description="The relative copy class for HGVS duplications and deletions represented as Relative Copy Number Variation."),  # noqa: E501
    untranslatable_returns_text: bool = Query(False, description=untranslatable_descr)
) -> NormalizeService:
    """Normalize and translate a human readable variation description to a single
        VRSATILE Variation Descriptor. Performs fully-justified allele normalization.
        Will liftover to GRCh38 and aligns to a priority transcript. Will make
        inferences about the query.

    :param str q: Human readable variation description on GRCh37 or GRCh38 assembly
    :param Optional[HGVSDupDelModeEnum] hgvs_dup_del_mode:
        Must be: `default`, `absolute_cnv`, `relative_cnv`, `repeated_seq_expr`,
        `literal_seq_expr`. This parameter determines how to interpret HGVS dup/del
        expressions in VRS.
    :param Optional[int] baseline_copies: Baseline copies for HGVS duplications and
        deletions. Required when `hgvs_dup_del_mode` is set to `absolute_cnv`.
    :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        for HGVS duplications and deletions represented as Relative Copy Number
        Variation. If not set, will use default `relative_copy_class` for query.
    :param bool untranslatable_returns_text: `True` return VRS Text Object when
        unable to translate or normalize query. `False` return `None` when
        unable to translate or normalize query.
    :return: NormalizeService for variation
    """
    normalize_resp = await query_handler.normalize_handler.normalize(
        unquote(q), hgvs_dup_del_mode=hgvs_dup_del_mode,
        baseline_copies=baseline_copies, relative_copy_class=relative_copy_class,
        untranslatable_returns_text=untranslatable_returns_text)
    return normalize_resp


@app.get("/variation/translate_identifier",
         summary="Given an identifier, use SeqRepo to return a list of aliases.",
         response_description="A response to a validly-formed query.",
         response_model=TranslateIdentifierService,
         description="Return list of aliases for an identifier",
         tags=[Tags.SEQREPO]
         )
def translate_identifier(
        identifier: str = Query(..., description="The identifier to find aliases for"),
        target_namespaces: Optional[str] = Query(None, description="The namespaces of the aliases, separated by commas")  # noqa: E501
) -> TranslateIdentifierService:
    """Return data containing identifier aliases.

    :param str identifier: The identifier to find aliases for
    :param Optional[str] target_namespaces: The namespaces of the aliases,
        separated by commas
    :return: TranslateIdentifierService data
    """
    aliases = []
    warnings = None
    identifier = identifier.strip()
    try:
        aliases = query_handler._seqrepo_access.seqrepo_client.translate_identifier(
            identifier, target_namespaces=target_namespaces)
    except KeyError:
        warnings = [f"Identifier, {identifier}, does not exist in SeqRepo"]
    except Exception as e:
        warnings = [f"SeqRepo could not translate identifier, {identifier}:"
                    f" {e}"]

    return TranslateIdentifierService(
        identifier_query=identifier,
        warnings=warnings,
        aliases=aliases,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        ))


from_fmt_descr = "Format of input variation to translate. Must be one of `beacon`, " \
                 "`gnomad`, `hgvs`, or `spdi`"


@app.get("/variation/translate_from",
         summary="Given variation as beacon, gnomad, hgvs or spdi representation, "
                 "return VRS Allele object using vrs-python's translator class",
         response_description="A response to a validly-formed query.",
         description="Return VRS Allele object",
         response_model=TranslateFromService,
         tags=[Tags.VRS_PYTHON])
def vrs_python_translate_from(
    variation: str = Query(..., description="Variation to translate to VRS object."
                                            " Must be represented as either beacon, "
                                            "gnomad, hgvs, or spdi."),
    fmt: Optional[TranslateFromFormat] = Query(None, description=from_fmt_descr)
) -> TranslateFromService:
    """Given variation query, return VRS Allele object using vrs-python"s translator
        class

    :param str variation: Variation to translate to VRS object. Must be represented
        as either beacon, gnomad, hgvs, or spdi
    :param Optional[TranslateFromFormat] fmt: Format of variation. If not supplied,
        vrs-python will infer its format.
    :return: TranslateFromService containing VRS Allele object
    """
    variation_query = unquote(variation.strip())
    warnings = list()
    vrs_variation = None
    try:
        resp = query_handler._tlr.translate_from(variation_query, fmt)
    except (KeyError, ValueError, python_jsonschema_objects.validators.ValidationError) as e:  # noqa: E501
        warnings.append(f"vrs-python translator raised {type(e).__name__}: {e}")
    except HGVSError as e:
        warnings.append(f"hgvs raised {type(e).__name__}: {e}")
    except BioutilsError as e:
        warnings.append(f"bioutils raised {type(e).__name__}: {e}")
    else:
        vrs_variation = resp.as_dict()

    return TranslateFromService(
        query=TranslateFromQuery(variation=variation_query, fmt=fmt),
        warnings=warnings,
        variation=vrs_variation,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        ),
        vrs_python_meta_=VrsPythonMeta(
            version=pkg_resources.get_distribution("ga4gh.vrs").version
        )
    )


g_to_p_summary = "Given gnomad VCF, return VRSATILE object on protein coordinate."
g_to_p_response_description = "A response to a validly-formed query."
g_to_p_description = \
    "Return VRSATILE object on protein coordinate for variation provided."
q_description = "gnomad VCF to normalize to protein variation."


@app.get("/variation/gnomad_vcf_to_protein",
         summary=g_to_p_summary,
         response_description=g_to_p_response_description,
         description=g_to_p_description,
         response_model=NormalizeService,
         tags=[Tags.TO_PROTEIN_VARIATION]
         )
async def gnomad_vcf_to_protein(
    q: str = Query(..., description=q_description),
    untranslatable_returns_text: bool = Query(False, description=untranslatable_descr)
) -> NormalizeService:
    """Return Value Object Descriptor for variation on protein coordinate.

    :param str q: gnomad VCF to normalize to protein variation.
    :param bool untranslatable_returns_text: `True` return VRS Text Object when
        unable to translate or normalize query. `False` return `None` when
        unable to translate or normalize query.
    :return: NormalizeService for variation
    """
    q = unquote(q.strip())
    resp = await query_handler.gnomad_vcf_to_protein_handler.gnomad_vcf_to_protein(
        q, untranslatable_returns_text=untranslatable_returns_text)
    return resp


complement_descr = "This field indicates that a categorical variation is defined to " \
                   "include (false) or exclude (true) variation concepts matching " \
                   "the categorical variation."
hgvs_dup_del_mode_decsr = "This parameter determines how to interpret HGVS dup/del "\
                          "expressions in VRS. Must be one of: `default`, " \
                          "`absolute_cnv`, `relative_cnv`, `repeated_seq_expr`, " \
                          "`literal_seq_expr`."
relative_copy_class_descr = "The relative copy class. Only used when `fmt`=`hgvs` "\
                            "and Relative CNV."
baseline_copies_descr = "Baseline copies for duplication or deletion. Only used when "\
                        "`fmt`=`hgvs` and Absolute CNV.`"


@app.get("/variation/to_canonical_variation",
         summary="Given SPDI or HGVS, return VRSATILE Canonical Variation",
         response_description="A response to a validly-formed query.",
         description="Return VRSATILE Canonical Variation",
         response_model=ToCanonicalVariationService,
         tags=[Tags.TO_CANONICAL])
async def to_canonical_variation(
        q: str = Query(..., description="HGVS or SPDI query"),
        fmt: ToCanonicalVariationFmt = Query(...,
                                             description="Format of the input variation. Must be `spdi` or `hgvs`"),  # noqa: E501
        complement: bool = Query(False, description=complement_descr),
        do_liftover: bool = Query(False, description="Whether or not to liftover to "
                                  "GRCh38 assembly."),
        hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = Query(
            HGVSDupDelModeEnum.DEFAULT, description=hgvs_dup_del_mode_decsr),
        relative_copy_class: Optional[RelativeCopyClass] = Query(
            None, description=relative_copy_class_descr),
        baseline_copies: Optional[int] = Query(
            None, description=baseline_copies_descr),
        untranslatable_returns_text: bool = Query(
            False, description=untranslatable_descr)
) -> ToCanonicalVariationService:
    """Return categorical variation for canonical SPDI

    :param str q: HGVS or SPDI query
    :param ToCanonicalVariationFmt fmt: Format of the input variation. Must be
        `spdi` or `hgvs`.
    :param bool complement: This field indicates that a categorical variation
        is defined to include (false) or exclude (true) variation concepts matching the
        categorical variation. This is equivalent to a logical NOT operation on the
        categorical variation properties.
    :param Optional[HGVSDupDelModeEnum] hgvs_dup_del_mode: Determines how to interpret
        HGVS dup/del expressions in VRS. Must be one of: `default`, `absolute_cnv`,
        `relative_cnv`, `repeated_seq_expr`, `literal_seq_expr`
    :param Optional[RelativeCopyClass] relative_copy_class: Relative copy class.
        Only used when `fmt`=`hgvs` and Relative CNV.
    :param Optional[int] baseline_copies: Baseline copies number
        Only used when `fmt`=`hgvs` and Absolute CNV
    :param bool untranslatable_returns_text: `True` return VRS Text Object when
        unable to translate or normalize query. `False` return `None` when
        unable to translate or normalize query.
    :return: ToCanonicalVariationService for variation query
    """
    q = unquote(q)
    resp = await query_handler.to_canonical_handler.to_canonical_variation(
        q, fmt, complement, do_liftover=do_liftover,
        hgvs_dup_del_mode=hgvs_dup_del_mode, baseline_copies=baseline_copies,
        relative_copy_class=relative_copy_class,
        untranslatable_returns_text=untranslatable_returns_text)
    return resp


def _get_allele(request_body: Union[TranslateToQuery, TranslateToHGVSQuery],
                warnings: List) -> Optional[models.Allele]:
    """Return VRS allele object from request body. `warnings` will get updated if
    exceptions are raised

    :param Union[TranslateToQuery, TranslateToHGVSQuery] request_body: Request body
        containing `variation`
    :param List warnings: List of warnings
    :return: VRS Allele object if valid
    """
    allele = None
    try:
        allele = models.Allele(**request_body["variation"])
    except ValidationError as e:
        warnings.append(f"`allele` is not a valid VRS Allele: {e}")
    except python_jsonschema_objects.ValidationError as e:
        warnings.append(str(e))
    return allele


@app.post("/variation/translate_to",
          summary="Given VRS Allele object as a dict, return variation expressed as "
                  "queried format using vrs-python's translator class",
          response_description="A response to a validly-formed query.",
          description="Return variation in queried format representation. "
                      "Request body must contain `variation` and `fmt`. `variation` is"
                      " a VRS Allele object represented as a dict. `fmt` must be either"
                      " `spdi` or `hgvs`",
          response_model=TranslateToService,
          tags=[Tags.VRS_PYTHON])
async def vrs_python_translate_to(
        request_body: TranslateToQuery) -> TranslateToService:
    """Given VRS Allele object as a dict, return variation expressed as queried
    format using vrs-python's translator class

    :param TranslateToQuery request_body: Request body. `variation` is a VRS Allele
        object represented as a dict. `fmt` must be either `spdi` or `hgvs`
    :return: TranslateToService containing variation represented as fmt representation
        if valid VRS Allele, and warnings if found
    """
    query = request_body
    request_body = request_body.dict(by_alias=True)
    warnings = list()

    allele = _get_allele(request_body, warnings)

    variations = list()
    if allele:
        try:
            variations = query_handler._tlr.translate_to(allele, request_body["fmt"])
        except ValueError as e:
            warnings.append(f"vrs-python translator raised {type(e).__name__}: {e}")

    return TranslateToService(
        query=query,
        warnings=warnings,
        variations=variations,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        ),
        vrs_python_meta_=VrsPythonMeta(
            version=pkg_resources.get_distribution("ga4gh.vrs").version
        )
    )


to_hgvs_descr = "Return variation as HGVS expressions. Request body must"\
                " contain `variation`, a VRS Allele object represented as a dict. "\
                "Can include optional parameter `namespace`. If `namespace` is not"\
                " None, returns HGVS strings for the specified namespace. If "\
                "`namespace` is None, returns HGVS strings for all alias translations."


@app.post("/variation/vrs_allele_to_hgvs",
          summary="Given VRS Allele object as a dict, return HGVS expression(s)",
          response_description="A response to a validly-formed query.",
          description=to_hgvs_descr,
          response_model=TranslateToService,
          tags=[Tags.VRS_PYTHON])
async def vrs_python_to_hgvs(
        request_body: TranslateToHGVSQuery) -> TranslateToService:
    """Given VRS Allele object as a dict, return variation expressed as HGVS
        expression(s)

    :param TranslateToHGVSQuery request_body: Request body. `variation` is a VRS Allele
        object represented as a dict. Can provide optional parameter `namespace`.
        If `namespace` is not None, returns HGVS strings for the specified namespace.
        If `namespace` is None, returns HGVS strings for all alias translations.
    :return: TranslateToService containing variation represented as HGVS representation
        if valid VRS Allele, and warnings if found
    """
    query = request_body
    request_body = request_body.dict(by_alias=True)
    warnings = list()

    allele = _get_allele(request_body, warnings)

    variations = list()
    if allele:
        try:
            variations = query_handler._tlr._to_hgvs(
                allele, namespace=request_body.get("namespace") or "refseq")
        except ValueError as e:
            warnings.append(f"vrs-python translator raised {type(e).__name__}: {e}")

    return TranslateToService(
        query=query,
        warnings=warnings,
        variations=variations,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        ),
        vrs_python_meta_=VrsPythonMeta(
            version=pkg_resources.get_distribution("ga4gh.vrs").version
        )
    )


@app.get("/variation/hgvs_to_absolute_copy_number",
         summary="Given HGVS expression, return VRS Absolute Copy Number Variation",
         response_description="A response to a validly-formed query.",
         description="Return VRS Absoluter Copy Number Variation",
         response_model=HgvsToAbsoluteCopyNumberService,
         tags=[Tags.TO_COPY_NUMBER_VARIATION])
async def hgvs_to_absolute_copy_number(
    hgvs_expr: str = Query(..., description="Variation query"),
    baseline_copies: Optional[int] = Query(
        None, description="Baseline copies for duplication"),
    do_liftover: bool = Query(False, description="Whether or not to liftover "
                              "to GRCh38 assembly."),
    untranslatable_returns_text: bool = Query(False, description=untranslatable_descr)
) -> HgvsToAbsoluteCopyNumberService:
    """Given hgvs expression, return absolute copy number variation

    :param str hgvs_expr: HGVS expression
    :param Optional[int] baseline_copies: Baseline copies number
    :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
    :param bool untranslatable_returns_text: `True` return VRS Text Object when
        unable to translate or normalize query. `False` return `None` when
        unable to translate or normalize query.
    :return: HgvsToAbsoluteCopyNumberService
    """
    resp = await query_handler.to_copy_number_handler.hgvs_to_absolute_copy_number(
        unquote(hgvs_expr.strip()), baseline_copies, do_liftover,
        untranslatable_returns_text=untranslatable_returns_text)
    return resp


@app.get("/variation/hgvs_to_relative_copy_number",
         summary="Given HGVS expression, return VRS Relative Copy Number Variation",
         response_description="A response to a validly-formed query.",
         description="Return VRS Relative Copy Number Variation",
         response_model=HgvsToRelativeCopyNumberService,
         tags=[Tags.TO_COPY_NUMBER_VARIATION])
async def hgvs_to_relative_copy_number(
    hgvs_expr: str = Query(..., description="Variation query"),
    relative_copy_class: RelativeCopyClass = Query(
        ..., description="The relative copy class"),
    do_liftover: bool = Query(False, description="Whether or not to liftover "
                              "to GRCh38 assembly."),
    untranslatable_returns_text: bool = Query(False, description=untranslatable_descr)
) -> HgvsToRelativeCopyNumberService:
    """Given hgvs expression, return relative copy number variation

    :param str hgvs_expr: HGVS expression
    :param RelativeCopyClass relative_copy_class: Relative copy class
    :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
    :param bool untranslatable_returns_text: `True` return VRS Text Object when
        unable to translate or normalize query. `False` return `None` when
        unable to translate or normalize query.
    :return: HgvsToRelativeCopyNumberService
    """
    resp = await query_handler.to_copy_number_handler.hgvs_to_relative_copy_number(
        unquote(hgvs_expr.strip()), relative_copy_class, do_liftover,
        untranslatable_returns_text=untranslatable_returns_text)
    return resp


assembly_descr = "Assembly. If `accession` is set, will ignore `assembly` and `chr`. "\
                 "If `accession` not set, must provide both `assembly` and `chr`."
chr_descr = "Chromosome. Must set when `assembly` is set."
accession_descr = "Accession. If `accession` is set, will ignore `assembly` and "\
                  "`chr`. If `accession` not set, must provide both `assembly` and `chr`."  # noqa: E501
start_descr = "Start position as residue coordinate"
end_descr = "End position as residue coordinate"
total_copies_descr = "Total copies for Absolute Copy Number variation object"


@app.get("/variation/parsed_to_abs_cnv",
         summary="Given parsed ClinVar Copy Number Gain/Loss components, return "
         "VRS Absolute Copy Number Variation",
         response_description="A response to a validly-formed query.",
         description="Return VRS Absolute Copy Number Variation",
         response_model=ParsedToAbsCnvService,
         tags=[Tags.TO_COPY_NUMBER_VARIATION]
         )
def parsed_to_abs_cnv(
    assembly: Optional[ClinVarAssembly] = Query(None, description=assembly_descr),
    chr: Optional[str] = Query(None, description=chr_descr),
    accession: Optional[str] = Query(None, description=accession_descr),
    start: int = Query(..., description=start_descr),
    end: int = Query(..., description=end_descr),
    total_copies: int = Query(..., description=total_copies_descr),
    untranslatable_returns_text: bool = Query(False, description=untranslatable_descr)
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
    :return: Tuple containing Absolute Copy Number variation and list of warnings
    """
    resp = query_handler.to_copy_number_handler.parsed_to_abs_cnv(
        start, end, total_copies, assembly, chr, accession,
        untranslatable_returns_text=untranslatable_returns_text)
    return resp
