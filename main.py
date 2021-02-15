"""Main application for FastAPI."""
from fastapi import FastAPI, Query
from fastapi.openapi.utils import get_openapi
from varlexapp.classifiers import Classify
from varlexapp.tokenizers import Tokenize
from varlexapp.validators import Validate
from varlexapp.translators import Translate
from varlexapp.data_sources import SeqRepoAccess, TranscriptMappings
from varlexapp.models import TranslationResponse, ValidationResponse
from varlexapp.schemas import TranslationResponseSchema, \
    ValidationResponseSchema, ClassificationResponseSchema, TokenResponseSchema

app = FastAPI(docs_url='/variant', openapi_url='/variant/openapi.json')
tokenizer = Tokenize('varlexapp/data/gene_symbols.txt')
classifier = Classify()
seq_repo_access = SeqRepoAccess('varlexapp/data/seqrepo/latest')
transcript_mappings = \
    TranscriptMappings('varlexapp/data/transcript_mapping.tsv')
validator = Validate(seq_repo_access, transcript_mappings)
translator = Translate(seq_repo_access)


def custom_openapi():
    """Generate custom fields for OpenAPI response."""
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="The VICC Variant Normalizer",
        version="0.1.0",
        description="VarLex - 'variant lexicon' prototype for "
                    "variant normalization.",
        routes=app.routes
    )

    openapi_schema['info']['contact'] = {
        "name": "Alex H. Wagner",
        "email": "Alex.Wagner@nationwidechildrens.org"
    }
    app.openapi_schema = openapi_schema
    return app.openapi_schema


app.openapi = custom_openapi


token_summary = "Tokenize a given query."
token_description = "Return the detected tokens for a given query."
token_response_description = "A  response to a validly-formed query."
q_description = "Query to tokenize."


@app.get('/variant/tokenization',
         summary=token_summary,
         response_description=token_response_description,
         response_model=TokenResponseSchema,
         description=token_description)
def tokenize(q: str = Query(..., description="Tokenize a given variant.")):
    """Return the tokens that VarLex detected in the query.

    :param str q: The query to search on
    """
    tokens = tokenizer.perform(q)
    return TokenResponseSchema(search_term=q, tokens=tokens)


classify_summary = "Classify a given query."
classify_description = ("Return the detected classifications with associated "
                        "confidence levels for a given variant.")
classify_response_description = "A  response to a validly-formed query."
q_description = "Variant to classify."


@app.get('/variant/classification',
         summary=classify_summary,
         response_description=classify_response_description,
         response_model=ClassificationResponseSchema,
         description=classify_description)
def classify(q: str = Query(..., description="Classify a given variant.")):
    """Return the classifications that VarLex detected in the query.

    :param str q: The query to search on
    """
    tokens = tokenizer.perform(q)
    classifications = classifier.perform(tokens)
    return ClassificationResponseSchema(search_term=q,
                                        classifications=classifications)


validate_summary = "Validate a given variant."
validate_description = ("Return the determined validation status "
                        "for a given variant.")
validate_response_description = "A  response to a validly-formed query."
q_description = "Variant to validate."


@app.get('/variant/validation',
         summary=validate_summary,
         response_description=validate_response_description,
         # response_model=ValidationResponseSchema,
         description=validate_description)
def validate(q: str = Query(..., description="Validate a given variant.")):
    """Return the validation status that VarLex determined for a given variant.

    :param str q: The variant to search on
    """
    tokens = tokenizer.perform(q)
    classifications = classifier.perform(tokens)
    res = validator.perform(classifications)
    resp = ValidationResponse(q, validation_summary=res)
    return ValidationResponseSchema().dump(resp)


translate_summary = "Translate a given variant."
translate_description = ("Translate a variant into VR-Spec (VRS) compatible "
                         "representations.")
translate_response_description = "A  response to a validly-formed query."
q_description = "Variant to translate."


@app.get('/variant/translate',
         summary=translate_summary,
         response_description=translate_response_description,
         # response_model=TranslationResponseSchema,
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

    resp = TranslationResponse(q, translations)

    return TranslationResponseSchema().dump(resp)
