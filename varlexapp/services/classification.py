from flask import (Blueprint, request)


from ..models import ClassificationResponse
from ..schemas import ClassificationResponseSchema
from ..classifiers import Classify
from ..tokenizers import Tokenize

bp = Blueprint('classify', __name__, url_prefix='/classifications')

tokenizer = Tokenize('varlexapp/data/gene_symbols.txt')
classifier = Classify()

@bp.route('/')
def get_classifications():
    """Returns the classifications that varlex detected in your query
    ---
    get:
        description: Get the detected classifications with associated confidence levels for a variant
        parameters:
            - name: q
              in: query
              required: true
              description: variant to classify
              schema:
                type: string
        responses:
            200:
                description: OK
                content:
                    application/json:
                        schema:  ClassificationResponseSchema
    """

    search_term = request.args.get('q')
    tokens = tokenizer.perform(search_term)
    classifications = classifier.perform(tokens)
    resp = ClassificationResponse(search_term=search_term, classifications=classifications)
    return ClassificationResponseSchema().dump(resp)
