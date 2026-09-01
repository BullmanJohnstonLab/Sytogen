import os
import logging

from flask import Flask, jsonify
from werkzeug.exceptions import RequestEntityTooLarge
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Hard cap on request body size (applies to every route, including file
# uploads). Without this, Flask/Werkzeug will happily buffer a request of
# any size before a view function ever runs - so an unbounded upload was
# a disk/memory exhaustion vector even on endpoints that later validate
# the *parsed* content (e.g. the 8 Mb construct-size check in api.py
# only runs after the file has already been fully saved to disk).
#
# 25 MB comfortably covers a full bacterial genome in GenBank format
# (a few MB, even with dense annotations) with headroom to spare, while
# still bounding the worst case. Override with SYTOGEN_MAX_UPLOAD_BYTES
# for deployments that need a different limit.
MAX_UPLOAD_BYTES = int(os.environ.get("SYTOGEN_MAX_UPLOAD_BYTES", 25 * 1024 * 1024))

# Rate limiting: 50 requests per hour per IP for heavy compute endpoints,
# 1000 per hour for light endpoints. Override with environment variables.
HEAVY_RATE_LIMIT = os.environ.get("SYTOGEN_HEAVY_RATE_LIMIT", "50 per hour")
LIGHT_RATE_LIMIT = os.environ.get("SYTOGEN_LIGHT_RATE_LIMIT", "1000 per hour")


def create_app():
    app = Flask(__name__)
    app.config["MAX_CONTENT_LENGTH"] = MAX_UPLOAD_BYTES
    
    # Initialize rate limiter
    limiter = Limiter(
        app=app,
        key_func=get_remote_address,
        default_limits=[LIGHT_RATE_LIMIT],
        storage_uri="memory://"
    )

    @app.errorhandler(RequestEntityTooLarge)
    def _handle_upload_too_large(e):
        limit_mb = MAX_UPLOAD_BYTES / (1024 * 1024)
        logger.warning(f"Upload rejected: exceeded {limit_mb:.0f} MB limit from {get_remote_address()}")
        return jsonify(
            error=f"Upload too large. The maximum request size is {limit_mb:.0f} MB."
        ), 413

    # Register web (HTML) routes
    from .web import web
    app.register_blueprint(web)

    # Register API routes
    from .api import api, start_job_sweeper
    app.register_blueprint(api, url_prefix="/api")
    
    # Store limiter reference for api.py to use
    app.limiter = limiter
    
    start_job_sweeper()
    logger.info("SyToGen app initialized")

    return app