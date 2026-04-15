# Copilot instructions for SyToGen

- Root app entrypoint: `run.py`. `python run.py` starts the Flask server in development mode.
- App factory is in `sytogen/__init__.py`. It registers two blueprints:
  - `sytogen/web.py` for HTML pages and form routes
  - `sytogen/api.py` for backend API routes under `/api`

- Primary feature implementation lives in `sytogen/scripts/motiffinder_backend.py`.
  - Route logic should remain in `sytogen/api.py`
  - Pure motif parsing, GFF3 loading, circular-sequence handling, and TSV/GFF3 output generation belong in the backend module

- Front-end behavior is driven by `sytogen/templates/motiffinder.html`.
  - It sends `POST /api/motiffinder/run` with multipart form data
  - Required fields: `sequence_file`, `motif_file`, `source_type`
  - Optional field: `annotation_file` only for FASTA source files
  - The JS client uses format sniffing to detect `genbank`, `fasta`, or `gff3`

- Request conventions to preserve:
  - `source_type` must be `genbank` or `fasta`
  - `annotation_file` is only included when `source_type == "fasta"`
  - GenBank uploads may include feature annotations directly in the file

- Dependencies are declared in `requirements.txt`; key libraries:
  - `Flask==3.1.2`
  - `biopython>=1.80`
  - `pytest>=7.0` (no test modules currently exist in the repo root)

- If adding new pages or UI routes, update `sytogen/web.py` and the corresponding Jinja templates under `sytogen/templates/`.
- If adding new API behavior, keep the public path prefix as `/api` and make changes in `sytogen/api.py`.
- There is no external database or cloud service in the current repo; the app is a self-contained Flask web service.

- Debugging note: `sytogen/api.py` uses `abort(400, ...)` on invalid uploads, and `motiffinder.html` displays the returned plain-text error message.
