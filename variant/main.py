"""Main application for FastAPI."""
from fastapi import FastAPI, Query
from fastapi.openapi.utils import get_openapi
from variant.to_vrs import ToVRS
from variant.schemas import ToVRSService, NormalizeService, ServiceMeta
from variant.normalize import Normalize
from .version import __version__
from datetime import datetime
import html

app = FastAPI(docs_url='/variant', openapi_url='/variant/openapi.json')

to_vrs = ToVRS()
normalizer = Normalize()


def custom_openapi():
    """Generate custom fields for OpenAPI response."""
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="The VICC Variant Normalizer",
        version=__version__,
        description="Services and guidelines for normalizing variants.",
        routes=app.routes
    )

    openapi_schema['info']['contact'] = {
        "name": "Alex H. Wagner",
        "email": "Alex.Wagner@nationwidechildrens.org"
    }
    app.openapi_schema = openapi_schema
    return app.openapi_schema


app.openapi = custom_openapi

translate_summary = "Translate a variant to a VRS compatible object."
translate_description = ("Translate a variant into Variation Representation "
                         "Specification (VRS) compatible representations.")
translate_response_description = "A  response to a validly-formed query."
q_description = "Variant to translate."


@app.get('/variant/toVRS',
         summary=translate_summary,
         response_description=translate_response_description,
         response_model=ToVRSService,
         description=translate_description)
def translate(q: str = Query(..., description=q_description)):
    """Return a VRS-like representation of all validated variants for the search term.  # noqa: E501, D400

    :param str q: The variant to search on
    :return: TranslationResponseSchema for variant
    """
    validations, warnings = to_vrs.get_validations(html.unescape(q))
    translations, warnings = to_vrs.get_translations(validations, warnings)

    return ToVRSService(
        search_term=q,
        variants=translations,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        ),
        warnings=warnings
    )


normalize_summary = 'Given variant, return Value Object Descriptor.'
normalize_response_description = 'A response to a validly-formed query.'
normalize_description = 'Return Value Object Descriptor for variant provided.'
q_description = 'Variant to normalize.'


@app.get('/variant/normalize',
         summary=normalize_summary,
         response_description=normalize_response_description,
         response_model=NormalizeService,
         description=normalize_description,
         response_model_exclude_none=True
         )
def normalize(q: str = Query(..., description=q_description)):
    """Return Value Object Descriptor for variant.

    :param q: Variant to normalize
    :return: NormalizeService for variant
    """
    validations, warnings = to_vrs.get_validations(q)
    normalize_resp = normalizer.normalize(html.unescape(q),
                                          validations,
                                          to_vrs.amino_acid_cache,
                                          warnings)

    return NormalizeService(
        variant_query=q,
        variation_descriptor=normalize_resp,
        warnings=normalizer.warnings if normalizer.warnings else None,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        )
    )
