from flask import (Blueprint, request)

from ..data_sources import SeqRepoAccess
from ..classifiers import Classify
from ..tokenizers import Tokenize
from ..validators import Validate

from ..models import ValidationResponse
from ..schemas import ValidationResponseSchema

from ..data_sources import SeqRepoAccess

bp = Blueprint('validate', __name__, url_prefix='/validations')

tokenizer = Tokenize('varlexapp/data/gene_symbols.txt')
classifier = Classify()
seq_repo_access = SeqRepoAccess('varlexapp/data/seqrepo/latest')
validator = Validate(seq_repo_access)

@bp.route('/')
def get_validation_status():
    """Returns the validation status that varlex determined for your variant
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
                        schema:  ValidationResponseSchema
    """

    search_term = request.args.get('q')
    tokens = tokenizer.perform(search_term)
    classifications = classifier.perform(tokens)
    res = validator.perform(classifications)

    resp = ValidationResponse(search_term=search_term, validation_summary=res)
    return ValidationResponseSchema().dump(resp)
