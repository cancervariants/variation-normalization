from flask import (Blueprint, request)

bp = Blueprint('tokenize', __name__, url_prefix='/tokens')

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
    return { "tokens": [ request.args.get("q") ] }
