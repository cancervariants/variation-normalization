from flask import (Blueprint, request)

from ..data_sources import SeqRepoAccess, TranscriptMappings
from ..classifiers import Classify
from ..tokenizers import Tokenize
from ..validators import Validate
from ..translators import Translate

from ..models import TranslationResponse
from ..schemas import TranslationResponseSchema

from ..data_sources import SeqRepoAccess

bp = Blueprint('translate', __name__, url_prefix='/translations')

tokenizer = Tokenize('varlexapp/data/gene_symbols.txt')
classifier = Classify()
seq_repo_access = SeqRepoAccess('varlexapp/data/seqrepo/latest')
transcript_mappings = TranscriptMappings('varlexapp/data/transcript_mapping.tsv')
validator = Validate(seq_repo_access, transcript_mappings)
translator = Translate(seq_repo_access)

@bp.route('/')
def get_translations():
    """Returns returns a VRS-like representation of all validated variants for the search term
    ---
    get:
        description: Classify a variant (q) and provide validation status for it
        parameters:
            - name: q
              in: query
              required: true
              description: variant to validate
              schema:
                type: string
        responses:
            200:
                description: OK
                content:
                    application/json:
                        schema:  TranslationResponseSchema
    """

    search_term = request.args.get('q')
    tokens = tokenizer.perform(search_term)
    classifications = classifier.perform(tokens)
    validations = validator.perform(classifications)

    translations = []
    for valid_variant in validations.valid_results:
        translations.append(translator.perform(valid_variant))

    resp = TranslationResponse(search_term, translations)
    return TranslationResponseSchema().dump(resp)
