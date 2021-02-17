"""Main application for FastAPI."""
from fastapi import FastAPI, Query
from fastapi.openapi.utils import get_openapi
from variant.classifiers import Classify
from variant.tokenizers import Tokenize
from variant.validators import Validate
from variant.translators import Translate
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant.schemas import TranslationResponseSchema

app = FastAPI(docs_url='/variant', openapi_url='/variant/openapi.json')
tokenizer = Tokenize('variant/data/gene_symbols.txt')
classifier = Classify()
seq_repo_access = SeqRepoAccess('variant/data/seqrepo/latest')
transcript_mappings = \
    TranscriptMappings('variant/data/transcript_mapping.tsv')
validator = Validate(seq_repo_access, transcript_mappings)
translator = Translate(seq_repo_access)


def custom_openapi():
    """Generate custom fields for OpenAPI response."""
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="The VICC Variant Normalizer",
        version="0.1.0",
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
translate_description = ("Translate a variant into VR-Spec (VRS) compatible "
                         "representations.")
translate_response_description = "A  response to a validly-formed query."
q_description = "Variant to translate."


@app.get('/variant/toVRS',
         summary=translate_summary,
         response_description=translate_response_description,
         response_model=TranslationResponseSchema,
         description=translate_description)
def translate(q: str = Query(..., description=q_description)):
    """Return a VRS-like representation of all validated variants for the search term.  # noqa: E501, D400

    :param str q: The variant to search on
    """
    tokens = tokenizer.perform(q)
    classifications = classifier.perform(tokens)
    validations = validator.perform(classifications)

    translations = []
    for valid_variant in validations.valid_results:
        translations.append(translator.perform(valid_variant))
    return TranslationResponseSchema(
        search_term=q,
        variants=translations
    )
