"""Main application for FastAPI."""
from typing import Optional
import pkg_resources
from fastapi import FastAPI, Query
from fastapi.openapi.utils import get_openapi
import python_jsonschema_objects
from variation.schemas import ToVRSService, NormalizeService, ServiceMeta
from .schemas.vrs_python_translator_schema import TranslateFromFormat, \
    TranslateFromService, TranslateFromQuery, VrsPythonMeta
from .version import __version__
from datetime import datetime
import html
from variation.query import QueryHandler
from variation.schemas.normalize_response_schema \
    import HGVSDupDelMode as HGVSDupDelModeEnum, TranslateIdentifierService, \
    CanonicalSPDIToCategoricalVariationService
from hgvs.exceptions import HGVSError
from bioutils.exceptions import BioutilsError

app = FastAPI(
    docs_url="/variation",
    openapi_url="/variation/openapi.json",
    swagger_ui_parameters={"tryItOutEnabled": True}
)
query_handler = QueryHandler()


def custom_openapi():
    """Generate custom fields for OpenAPI response."""
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="The VICC Variation Normalizer",
        version=__version__,
        description="Services and guidelines for normalizing variations.",
        routes=app.routes
    )

    openapi_schema['info']['contact'] = {
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


@app.get('/variation/toVRS',
         summary=translate_summary,
         response_description=translate_response_description,
         response_model=ToVRSService,
         description=translate_description,
         response_model_exclude_none=True
         )
def to_vrs(q: str = Query(..., description=q_description)):
    """Return a VRS-like representation of all validated variations for the search term.  # noqa: E501, D400

    :param str q: The variation to translate
    :return: ToVRSService model for variation
    """
    translations, warnings = query_handler.to_vrs(html.unescape(q))

    return ToVRSService(
        search_term=q,
        variations=translations,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        ),
        warnings=warnings
    )


normalize_summary = 'Given variation, return VRSATILE compatible object.'
normalize_response_description = 'A response to a validly-formed query.'
normalize_description = \
    'Return VRSATILE compatible object for variation provided.'
q_description = 'Variation to normalize.'
hgvs_dup_del_mode_decsr = ('Must be one of: `default`, `cnv`, '
                           '`repeated_seq_expr`, `literal_seq_expr`. This'
                           ' parameter determines how to interpret HGVS '
                           'dup/del expressions in VRS.')


@app.get('/variation/normalize',
         summary=normalize_summary,
         response_description=normalize_response_description,
         response_model=NormalizeService,
         description=normalize_description,
         response_model_exclude_none=True
         )
def normalize(q: str = Query(..., description=q_description),
              hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = Query(
                  None, description=hgvs_dup_del_mode_decsr)):
    """Return Value Object Descriptor for variation.

    :param str q: Variation to normalize
    :param Optional[HGVSDupDelModeEnum] hgvs_dup_del_mode:
        Must be: `default`, `cnv`, `repeated_seq_expr`, `literal_seq_expr`.
        This parameter determines how to interpret HGVS dup/del expressions
        in VRS.
    :return: NormalizeService for variation
    """
    normalize_resp = query_handler.normalize(
        html.unescape(q), hgvs_dup_del_mode=hgvs_dup_del_mode)
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


@app.get('/variation/translate_identifier',
         summary='Given an identifier, use SeqRepo to return a list of aliases.',  # noqa: E501
         response_description='A response to a validly-formed query.',
         response_model=TranslateIdentifierService,
         description='Return list of aliases for an identifier'
         )
def translate_identifier(
        identifier: str = Query(..., description='The identifier to find aliases for'),  # noqa: E501
        target_namespaces: Optional[str] = Query(None, description='The namespaces of the aliases, separated by commas')  # noqa: E501
) -> TranslateIdentifierService:
    """Return data containing identifier aliases.

    :param str identifier: The identifier to find aliases for
    :param Optional[str] target_namespaces: The namespaces of the aliases,
        separated by commas
    :return: TranslateIdentifierService data
    """
    aliases = []
    warnings = None
    try:
        aliases = query_handler.seqrepo_access.seq_repo_client.translate_identifier(  # noqa: E501
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


g_to_p_summary = 'Given gnomad VCF, return VRSATILE compatible object on ' \
                 'protein coordinate.'
g_to_p_response_description = 'A response to a validly-formed query.'
g_to_p_description = \
    'Return VRSATILE compatible object on protein coordinate for ' \
    'variation provided.'
q_description = 'gnomad VCF to normalize to protein variation.'


@app.get('/variation/gnomad_vcf_to_protein',
         summary=g_to_p_summary,
         response_description=g_to_p_response_description,
         description=g_to_p_description,
         response_model=NormalizeService,
         response_model_exclude_none=True
         )
def gnomad_vcf_to_protein(q: str = Query(..., description=q_description)):
    """Return Value Object Descriptor for variation on protein coordinate.

    :param str q: gnomad VCF to normalize to protein variation.
    :return: NormalizeService for variation
    """
    q = html.unescape(q.strip())
    resp, warnings = query_handler.gnomad_vcf_to_protein(q)

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


@app.get("/variation/canonical_spdi_to_categorical_variation",
         summary="Given canonical SPDI, return VRSATILE categorical variation object",
         response_description="A response to a validly-formed query.",
         description="Return VRSATILE categorical variation object",
         response_model=CanonicalSPDIToCategoricalVariationService,
         response_model_exclude_none=True)
def canonical_spdi_to_categorical_variation(
        q: str = Query(..., description="Canonical SPDI"),
        complement: bool = Query(False, description=complement_descr)
) -> CanonicalSPDIToCategoricalVariationService:
    """Return categorical variation for canonical SPDI

    :param str q: Canonical SPDI
    :param bool complement: This field indicates that a categorical variation
        is defined to include (false) or exclude (true) variation concepts matching the
        categorical variation. This is equivalent to a logical NOT operation on the
        categorical variation properties.
    :return: CanonicalSPDIToCategoricalVariationService for variation query
    """
    q = html.unescape(q.strip())
    resp, warnings = query_handler.canonical_spdi_to_categorical_variation(
        q, complement=complement)
    return CanonicalSPDIToCategoricalVariationService(
        canonical_spdi_query=q,
        categorical_variation=resp,
        warnings=warnings,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        )
    )


from_fmt_descr = "Format of input variation to translate. Must be one of `beacon`, " \
                 "`gnomad`, `hgvs`, or `spdi`"


@app.get("/variation/translate_from",
         summary="Given variation as beacon, gnomad, hgvs or spdi representation, "
                 "return VRS Allele object using vrs-python's translator class",
         response_description="A response to a validly-formed query.",
         description="Return VRS Allele object",
         response_model=TranslateFromService,
         response_model_exclude_none=True)
def vrs_python_translate_from(
    variation: str = Query(..., description="Variation to translate to VRS object."
                                            " Must be represented as either beacon, "
                                            "gnomad, hgvs, or spdi."),
    fmt: Optional[TranslateFromFormat] = Query(None, description=from_fmt_descr)
) -> TranslateFromService:
    """Given variation query, return VRS Allele object using vrs-python's translator
        class

    :param str variation: Variation to translate to VRS object. Must be represented
        as either beacon, gnomad, hgvs, or spdi
    :param Optional[TranslateFromFormat] fmt: Format of variation. If not supplied,
        vrs-python will infer its format.
    :return: TranslateFromService containing VRS Allele object
    """
    variation_query = html.unescape(variation)
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


# @app.post("/variation/translate_to",
#           summary="Given VRS Allele object as a dict, return variation expressed as "
#                   "queried format using vrs-python's translator class",
#           response_description="A response to a validly-formed query.",
#           description="Return variation in queried format representation. "
#                       "Request body must contain `allele` and `fmt`. `allele` is a "
#                       "VRS Allele object represented as a dict. `fmt` must be either "
#                       "`spdi` or `hgvs`",
#           response_model=TranslateToService,
#           response_model_exclude_none=True)
# async def vrs_python_translate_to(
#         request: Request) -> Union[ErrorResponse, TranslateToService]:
#     """Given VRS Allele object as a dict, return variation expressed as queried
#     format using vrs-python's translator class
#
#     :param Request request: Request body. Request body must contain `allele` and `fmt`.  # noqa
#         `allele` is a VRS Allele object represented as a dict. `fmt` must be either
#         `spdi` or `hgvs`
#     :return: ErrorResponse if invalid request body. Else, TranslateToService containing  # noqa
#         variation represented as fmt representation if valid VRS Allele
#     """
#     r = await request.json()
#     warnings = list()
#
#     allele_query = r.get("allele")
#     if not allele_query:
#         warnings.append("Missing `allele`. Must be VRS Allele represented as a dict")
#     else:
#         if not isinstance(allele_query, dict):
#             warnings.append("`allele` must be a dict")
#
#     if warnings:
#         return ErrorResponse(errors=warnings)
#
#     fmt_query = r.get("fmt")
#     if not fmt_query:
#         warnings.append("Missing `fmt`. Must be either `hgvs` or `spdi`")
#     else:
#         if not isinstance(fmt_query, str):
#             warnings.append("`fmt` must be a str")
#         else:
#             fmt_query = fmt_query.strip()
#             valid_fmts = [v.value for k, v in TranslateToFormat.__members__.items()]
#             if fmt_query not in valid_fmts:
#                 warnings.append(f"{fmt_query} not a valid fmt. "
#                                 f"Must be one of {valid_fmts}")
#
#     if warnings:
#         return ErrorResponse(errors=warnings)
#
#     allele = None
#     try:
#         allele = models.Allele(**r["allele"])
#     except ValidationError as e:
#         warnings.append(f"`allele` is not a valid VRS Allele: {e}")
#
#     variation = []
#     if allele:
#         try:
#             variation = query_handler.tlr.translate_to(allele, r["fmt"])
#         except ValueError as e:
#             warnings.append(f"vrs-python translator raised {type(e).__name__}: {e}")
#
#     return TranslateToService(
#         query=TranslateToQuery(variation=allele_query, fmt=fmt_query),
#         warnings=warnings,
#         variations=variation,
#         service_meta_=ServiceMeta(
#             version=__version__,
#             response_datetime=datetime.now()
#         ),
#         vrs_python_meta_=VrsPythonMeta(
#             version=pkg_resources.get_distribution("ga4gh.vrs").version
#         )
#     )
