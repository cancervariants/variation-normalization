from flask import (Blueprint, request)

from ..tokenizers import Tokenize

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
    """

    search_term = request.args.get('q')
    tokens = tokenizer.perform(search_term)

    return {
                'searchTerm': search_term,
                'tokens': [token.as_json() for token in tokens]
           }
