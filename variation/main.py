"""Main application for FastAPI."""
from enum import Enum
from typing import Dict, Optional
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
from variation.schemas.service_schema import ClinVarAssembly, ParsedToAbsCnvQuery, \
    ParsedToAbsCnvService
from .version import __version__
from .schemas.vrs_python_translator_schema import TranslateFromFormat, \
    TranslateFromService, TranslateFromQuery, TranslateToQuery, TranslateToService, \
    VrsPythonMeta


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

translate_summary = "Translate a variation to a VRS compatible object."
translate_description = ("Translate a variation into Variation Representation "
                         "Specification (VRS) compatible representations.")
translate_response_description = "A  response to a validly-formed query."
q_description = "Variation to translate."


@app.get("/variation/to_vrs",
         summary=translate_summary,
         response_description=translate_response_description,
         response_model=ToVRSService,
         description=translate_description,
         response_model_exclude_none=True
         )
async def to_vrs(q: str = Query(..., description=q_description)) -> ToVRSService:
    """Return a VRS-like representation of all validated variations for the search term.  # noqa: E501, D400

    :param str q: The variation to translate
    :return: ToVRSService model for variation
    """
    translations, warnings = await query_handler.to_vrs(unquote(q))

    return ToVRSService(
        search_term=q,
        variations=translations,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        ),
        warnings=warnings
    )


normalize_summary = "Given variation, return VRSATILE compatible object."
normalize_response_description = "A response to a validly-formed query."
normalize_description = \
    "Return VRSATILE compatible object for variation provided."
q_description = "Variation to normalize."
hgvs_dup_del_mode_decsr = ("Must be one of: `default`, `absolute_cnv`, `relative_cnv`, "
                           "`repeated_seq_expr`, `literal_seq_expr`. This"
                           " parameter determines how to interpret HGVS "
                           "dup/del expressions in VRS.")


@app.get("/variation/normalize",
         summary=normalize_summary,
         response_description=normalize_response_description,
         response_model=NormalizeService,
         description=normalize_description,
         response_model_exclude_none=True
         )
async def normalize(
    q: str = Query(..., description=q_description),
    hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = Query(
        HGVSDupDelModeEnum.DEFAULT, description=hgvs_dup_del_mode_decsr),
    baseline_copies: Optional[int] = Query(
        None, description="Baseline copies for HGVS duplications and deletions represented as Absolute Copy Number Variation"),  # noqa: E501
    relative_copy_class: Optional[RelativeCopyClass] = Query(
        None, description="The relative copy class for HGVS duplications and deletions represented as Relative Copy Number Variation.")  # noqa: E501
) -> NormalizeService:
    """Return Value Object Descriptor for variation.

    :param str q: Variation to normalize
    :param Optional[HGVSDupDelModeEnum] hgvs_dup_del_mode:
        Must be: `default`, `absolute_cnv`, `relative_cnv`, `repeated_seq_expr`,
        `literal_seq_expr`. This parameter determines how to interpret HGVS dup/del
        expressions in VRS.
    :param Optional[int] baseline_copies: Baseline copies for HGVS duplications and
        deletions
    :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        for HGVS duplications and deletions represented as Relative Copy Number
        Variation.
    :return: NormalizeService for variation
    """
    normalize_resp = await query_handler.normalize(
        unquote(q), hgvs_dup_del_mode=hgvs_dup_del_mode,
        baseline_copies=baseline_copies, relative_copy_class=relative_copy_class)
    warnings = query_handler.normalize_handler.warnings if \
        query_handler.normalize_handler.warnings else None

    return NormalizeService(
        variation_query=q,
        variation_descriptor=normalize_resp,
        warnings=warnings,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        )
    )


@app.get("/variation/translate_identifier",
         summary="Given an identifier, use SeqRepo to return a list of aliases.",  # noqa: E501
         response_description="A response to a validly-formed query.",
         response_model=TranslateIdentifierService,
         description="Return list of aliases for an identifier",
         tags=[Tags.SEQREPO]
         )
def translate_identifier(
        identifier: str = Query(..., description="The identifier to find aliases for"),  # noqa: E501
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
        aliases = query_handler.seqrepo_access.seqrepo_client.translate_identifier(  # noqa: E501
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
         response_model_exclude_none=True,
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
        resp = query_handler.tlr.translate_from(variation_query, fmt)
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


g_to_p_summary = "Given gnomad VCF, return VRSATILE compatible object on " \
                 "protein coordinate."
g_to_p_response_description = "A response to a validly-formed query."
g_to_p_description = \
    "Return VRSATILE compatible object on protein coordinate for " \
    "variation provided."
q_description = "gnomad VCF to normalize to protein variation."


@app.get("/variation/gnomad_vcf_to_protein",
         summary=g_to_p_summary,
         response_description=g_to_p_response_description,
         description=g_to_p_description,
         response_model=NormalizeService,
         response_model_exclude_none=True,
         tags=[Tags.TO_PROTEIN_VARIATION]
         )
async def gnomad_vcf_to_protein(
    q: str = Query(..., description=q_description)
) -> NormalizeService:
    """Return Value Object Descriptor for variation on protein coordinate.

    :param str q: gnomad VCF to normalize to protein variation.
    :return: NormalizeService for variation
    """
    q = unquote(q.strip())
    resp, warnings = await query_handler.gnomad_vcf_to_protein(q)

    return NormalizeService(
        variation_query=q,
        variation_descriptor=resp,
        warnings=warnings,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        )
    )


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
         summary="Given SPDI or HGVS, return VRSATILE canonical variation",
         response_description="A response to a validly-formed query.",
         description="Return VRSATILE canonical variation object",
         response_model=ToCanonicalVariationService,
         response_model_exclude_none=True,
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
            None, description=baseline_copies_descr)
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
    :return: ToCanonicalVariationService for variation query
    """
    q = unquote(q)
    variation, warnings = await query_handler.to_canonical_variation(
        q, fmt, complement, do_liftover=do_liftover,
        hgvs_dup_del_mode=hgvs_dup_del_mode, baseline_copies=baseline_copies,
        relative_copy_class=relative_copy_class)

    return ToCanonicalVariationService(
        query=q,
        canonical_variation=variation,
        warnings=warnings,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        )
    )


@app.post("/variation/translate_to",
          summary="Given VRS Allele object as a dict, return variation expressed as "
                  "queried format using vrs-python's translator class",
          response_description="A response to a validly-formed query.",
          description="Return variation in queried format representation. "
                      "Request body must contain `variation` and `fmt`. `variation` is"
                      " a VRS Allele object represented as a dict. `fmt` must be either"
                      " `spdi` or `hgvs`",
          response_model=TranslateToService,
          response_model_exclude_none=True,
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

    allele = None
    try:
        allele = models.Allele(**request_body["variation"])
    except ValidationError as e:
        warnings.append(f"`allele` is not a valid VRS Allele: {e}")
    except python_jsonschema_objects.ValidationError as e:
        warnings.append(str(e))

    variations = list()
    if allele:
        try:
            variations = query_handler.tlr.translate_to(allele, request_body["fmt"])
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
         summary="Given HGVS expression, return absolute copy number variation",
         response_description="A response to a validly-formed query.",
         description="Return VRS object",
         response_model=HgvsToAbsoluteCopyNumberService,
         response_model_exclude_none=True,
         tags=[Tags.TO_COPY_NUMBER_VARIATION])
async def hgvs_to_absolute_copy_number(
    hgvs_expr: str = Query(..., description="Variation query"),
    baseline_copies: Optional[int] = Query(
        None, description="Baseline copies for duplication"),
    do_liftover: bool = Query(False, description="Whether or not to liftover "
                              "to GRCh38 assembly.")
) -> HgvsToAbsoluteCopyNumberService:
    """Given hgvs expression, return absolute copy number variation

    :param str hgvs_expr: HGVS expression
    :param Optional[int] baseline_copies: Baseline copies number
    :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
    :return: HgvsToAbsoluteCopyNumberService
    """
    variation, warnings = await query_handler.hgvs_to_absolute_copy_number(
        unquote(hgvs_expr.strip()), baseline_copies, do_liftover)

    return HgvsToAbsoluteCopyNumberService(
        hgvs_expr=hgvs_expr,
        warnings=warnings,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        ),
        absolute_copy_number=variation
    )


@app.get("/variation/hgvs_to_relative_copy_number",
         summary="Given HGVS expression, return relative copy number variation",
         response_description="A response to a validly-formed query.",
         description="Return VRS object",
         response_model=HgvsToRelativeCopyNumberService,
         response_model_exclude_none=True,
         tags=[Tags.TO_COPY_NUMBER_VARIATION])
async def hgvs_to_relative_copy_number(
    hgvs_expr: str = Query(..., description="Variation query"),
    relative_copy_class: RelativeCopyClass = Query(
        ..., description="The relative copy class"),
    do_liftover: bool = Query(False, description="Whether or not to liftover "
                              "to GRCh38 assembly.")
) -> HgvsToRelativeCopyNumberService:
    """Given hgvs expression, return relative copy number variation

    :param str hgvs_expr: HGVS expression
    :param RelativeCopyClass relative_copy_class: Relative copy class
    :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
    :return: HgvsToRelativeCopyNumberService
    """
    variation, warnings = await query_handler.hgvs_to_relative_copy_number(
        unquote(hgvs_expr.strip()), relative_copy_class, do_liftover)

    return HgvsToRelativeCopyNumberService(
        hgvs_expr=hgvs_expr,
        warnings=warnings,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        ),
        relative_copy_number=variation
    )


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
         "absolute copy number variation",
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
    total_copies: int = Query(..., description=total_copies_descr)
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
    :return: Tuple containing Absolute Copy Number variation and list of warnings
    """
    variation, warnings = query_handler.parsed_to_abs_cnv(
        start, end, total_copies, assembly, chr, accession)
    query = ParsedToAbsCnvQuery(assembly=assembly, chr=chr, accession=accession,
                                start=start, end=end, total_copies=total_copies)
    return ParsedToAbsCnvService(
        query=query,
        absolute_copy_number=variation,
        warnings=warnings,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        )
    )
