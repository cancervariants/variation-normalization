import os

from flask import Flask
from flask_cors import CORS
from apispec import APISpec
from apispec_webframeworks.flask import FlaskPlugin
from apispec.ext.marshmallow import MarshmallowPlugin

def create_app(test_config=None):
    spec = APISpec(
            title="varlex",
            version="0.0.1",
            openapi_version="3.0.2",
            info=dict(description="Testing out the flask plugin"),
            plugins=[FlaskPlugin(), MarshmallowPlugin()]
    )

    app = Flask(__name__, instance_relative_config=True)
    CORS(app)
    app.config.from_mapping(SECRET_KEY='dev')

    if test_config is None:
        app.config.from_pyfile('config.py', silent=True)
    else:
        app.config.from_mapping(test_config)

    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    from .services import tokenization
    app.register_blueprint(tokenization.bp)

    with app.test_request_context():
        spec.path(view=tokenization.get_tokens)

    @app.route('/openapi.json')
    def get_openapi():
        return spec.to_dict()

    return app
