"""Main application for FastAPI."""
from typing import Optional
from fastapi import FastAPI, Query
from fastapi.openapi.utils import get_openapi
from variation.schemas import ToVRSService, NormalizeService, ServiceMeta
from .version import __version__
from datetime import datetime
import html
from variation.query import QueryHandler
from variation.schemas.normalize_response_schema \
    import HGVSDupDelMode as HGVSDupDelModeEnum, TranslateIdentifierService

app = FastAPI(docs_url='/variation', openapi_url='/variation/openapi.json')
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
        "email": "Alex.Wagner@nationwidechildrens.org"
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
