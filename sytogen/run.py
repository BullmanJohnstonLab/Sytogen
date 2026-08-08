import os

from sytogen import create_app

app = create_app()

if __name__ == "__main__":
    # Debug mode enables the Werkzeug interactive debugger, which allows
    # arbitrary code execution from the browser. It must stay off unless
    # explicitly requested for local development.
    debug = os.environ.get("FLASK_DEBUG", "0") == "1"
    app.run(debug=debug)
