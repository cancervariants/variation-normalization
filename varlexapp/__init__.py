import os

from flask import Flask
from apispec import APISpec
from apispec_webframeworks.flask import FlaskPlugin

def create_app(test_config=None):
    spec = APISpec(
            title="varlex",
            version="0.0.1",
            openapi_version="3.0.2",
            info=dict(description="Testing out the flask plugin"),
            plugins=[FlaskPlugin()]
    )

    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(SECRET_KEY='dev')

    if test_config is None:
        app.config.from_pyfile('config.py', silent=True)
    else:
        app.config.from_mapping(test_config)

    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    from . import tokenize
    app.register_blueprint(tokenize.bp)

    with app.test_request_context():
        spec.path(view=tokenize.get_tokens)

    @app.route('/swagger.json')
    def get_swagger():
        return spec.to_dict()

    return app
