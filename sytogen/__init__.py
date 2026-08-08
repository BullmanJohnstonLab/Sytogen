import os

from flask import Flask, jsonify
from werkzeug.exceptions import RequestEntityTooLarge

# Hard cap on request body size (applies to every route, including file
# uploads). Without this, Flask/Werkzeug will happily buffer a request of
# any size before a view function ever runs - so an unbounded upload was
# a disk/memory exhaustion vector even on endpoints that later validate
# the *parsed* content (e.g. the 20 kb construct-size check in api.py
# only runs after the file has already been fully saved to disk).
#
# 25 MB comfortably covers a full bacterial genome in GenBank format
# (a few MB, even with dense annotations) with headroom to spare, while
# still bounding the worst case. Override with SYTOGEN_MAX_UPLOAD_BYTES
# for deployments that need a different limit.
MAX_UPLOAD_BYTES = int(os.environ.get("SYTOGEN_MAX_UPLOAD_BYTES", 25 * 1024 * 1024))


def create_app():
    app = Flask(__name__)
    app.config["MAX_CONTENT_LENGTH"] = MAX_UPLOAD_BYTES

    @app.errorhandler(RequestEntityTooLarge)
    def _handle_upload_too_large(e):
        limit_mb = MAX_UPLOAD_BYTES / (1024 * 1024)
        return jsonify(
            error=f"Upload too large. The maximum request size is {limit_mb:.0f} MB."
        ), 413

    # Register web (HTML) routes
    from .web import web
    app.register_blueprint(web)

    # Register API routes
    from .api import api, start_job_sweeper
    app.register_blueprint(api, url_prefix="/api")
    start_job_sweeper()

    return app