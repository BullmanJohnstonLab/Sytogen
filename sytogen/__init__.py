from flask import Flask


def create_app():
    app = Flask(__name__)

    # Register web (HTML) routes
    from .web import web
    app.register_blueprint(web)

    # Register API routes
    from .api import api, start_job_sweeper
    app.register_blueprint(api, url_prefix="/api")
    start_job_sweeper()

    return app