# SyToGen

SyToGen (the **Sy**ngenicDNA **To**ol **Gen**erator) removes restriction–modification (RM)
target motifs from DNA constructs, so that a construct can be synthesized and transformed
into a bacterial host without being degraded by that host's own restriction–modification
system. It is host- and topology-aware: edits are chosen using strain-specific codon-usage
data, and both linear and circular constructs are supported.

SyToGen ships as a Flask web application with four linked tools, run in sequence:

| Tool | Page | What it does |
|---|---|---|
| **MyMotifs** | `/mymotif` | Parses REBASE/MIJAMP methylome output (or a manually curated list) into a standardized RM motif profile. |
| **CodonBias** | `/codon-bias` | Builds a strain-specific codon-usage table from an annotated genome. |
| **MotifFinder** | `/motiffinder` | Locates RM motif occurrences in a target construct, on both strands, accounting for circular topology. |
| **SyToGen** | `/sytogen` | The optimization engine: eliminates motifs via synonymous/neutral substitutions, ranked by motif destruction, avoidance of new motifs, host codon preference, and (optionally) GC preservation. Outputs the edited sequence, a full decision matrix, and a Gibson Assembly plan. |

A companion R package, [SytogenR](https://github.com/BullmanJohnstonLab/SytogenR), exposes
the same motif-parsing, codon-bias, motif-finding, and pipeline functions for use from R.

## Try it hosted

`[hosted app link]`

## Run it locally

```bash
git clone https://github.com/BullmanJohnstonLab/Sytogen.git
cd Sytogen
pip install -r requirements.txt
python run.py
```

This starts a local Flask dev server (`debug=True`) — by default at
`http://127.0.0.1:5000`. Open that address in a browser to reach the same pages listed
above.

**Requirements:** Python 3.9+ is recommended. Key dependencies (see `requirements.txt`
for exact pins): Flask, Biopython, pandas, NumPy, SciPy, scikit-learn, PuLP,
`dna_features_viewer`, `pydna`, Plotly.

## Using the API directly

Each tool is also reachable as a JSON API, which is what the web pages call under the
hood. All endpoints are under `/api`:

| Endpoint | Method | Purpose |
|---|---|---|
| `/api/mymotif/parse` | POST | Parse an uploaded motif file into a MyMotifs profile. |
| `/api/motiffinder/run` | POST | Run MotifFinder against a sequence + motif profile. |
| `/api/codonbias/run` | POST | Compute a codon-usage table from an uploaded genome. |
| `/api/sytogen/run` | POST | Run the SyToGen optimization synchronously. |
| `/api/sytogen/submit` | POST | Submit a SyToGen job asynchronously (long-running constructs). |
| `/api/status/<job_id>` | GET | Poll the status of a submitted job. |
| `/api/result/<job_id>` | GET | Retrieve the completed job's outputs. |

See the `/user-guide` and `/explained` pages on a running instance for input-format
details (accepted file types, required columns, etc.).

## Command-line tools

Install the project in an environment with the dependencies from
`requirements.txt`, then use the `sytogen` command:

```bash
pip install -e .
sytogen --version
sytogen mymotifs parse motifs.txt --output motifs.csv
sytogen codon-bias \
	--genome genome.gbk \
	--output-dir codon-bias-results
sytogen motif-finder \
	--sequence construct.gbk \
	--motifs motifs.csv \
	--topology circular \
	--output-dir motif-finder-results
sytogen run \
	--sequence construct.gbk \
	--codon-usage codon_usage.csv \
	--motifs motifs.csv \
	--topology circular \
	--output-dir sytogen-results \
	--zip sytogen-results.zip
```

`mymotifs parse` accepts CSV, TSV, REBASE, and MIJAMP motif files. `sytogen run`
accepts GenBank directly, or FASTA with `--source-type fasta --gff annotations.gff3`,
and writes the edited sequence, decision matrix, motif summaries, and validation
JSON to the output directory. Use `--assembly-plan` to include assembly fragments
and primers.

`codon-bias` also accepts `--fasta genome.fasta --gff annotations.gff3`. `motif-finder`
accepts GenBank or FASTA with optional GFF3 annotations, supports multiple sequence
records, and writes combined hit tables as TSV and GFF3 plus a per-record JSON summary.

## Project layout

```
Sytogen/
├── run.py                      # Flask entry point (python run.py)
├── requirements.txt
├── sytogen/
│   ├── __init__.py             # app factory (create_app)
│   ├── web.py                  # HTML page routes
│   ├── api.py                  # JSON API routes
│   ├── scripts/
│   │   ├── codon_bias_estimator.py
│   │   ├── motiffinder_backend.py
│   │   ├── sytogen_runner.py   # core optimization engine
│   │   └── visualization.py
│   └── templates/              # Jinja2 page templates
└── tests/
```

## Tests

```bash
pip install -r requirements.txt   # includes pytest
pytest
```

## Citation

If you use SyToGen, please cite: `[add citation once published]`.

## License

`[add license]`
