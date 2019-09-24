from flask import (Blueprint, request)

from ..tokenizers import Tokenize
from ..schemas import TokenResponseSchema
from ..models import TokenResponse

bp = Blueprint('tokenize', __name__, url_prefix='/tokens')

tokenizer = Tokenize('varlexapp/data/gene_symbols.txt')

@bp.route('/')
def get_tokens():
    """Returns the tokens that varlex detected in your query
    ---
    get:
        description: Get the detected tokens for a query "q"
        parameters:
            - name: q
              in: query
              required: true
              description: search term to tokenize
              schema:
                type: string
        responses:
            200:
                description: OK
                content:
                    application/json:
                        schema:  TokenResponseSchema
    """

    search_term = request.args.get('q')
    tokens = tokenizer.perform(search_term)

    resp = TokenResponse(search_term=search_term, tokens=tokens)

    return TokenResponseSchema().dump(resp)
