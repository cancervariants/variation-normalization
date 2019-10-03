from flask import (Blueprint, request)


from ..models import ClassificationResponse
from ..schemas import ClassificationResponseSchema

bp = Blueprint('classify', __name__, url_prefix='/classifications')


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
    resp = ClassificationResponse(search_term=search_term, classifications = [])
    return ClassificationResponseSchema().dump(resp)


