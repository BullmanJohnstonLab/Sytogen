import os
import tempfile
import warnings
from types import SimpleNamespace


class SyToGenPipelineError(RuntimeError):
    pass


def _first_existing(paths, *keys):
    for key in keys:

        path = paths.get(key)
        print("KEY:", key)
        print("PATH:", path)

        if path:
            print("EXISTS:", os.path.exists(path))
            
        if path and os.path.exists(path):
            return path
    return None

def count_motifs(seq, motifs):
    return sum(seq.count(m["motif"]) for m in motifs)


def _legacy_module():
    cache_dir = os.path.join(tempfile.gettempdir(), "sytogen-cache")
    matplotlib_dir = os.path.join(cache_dir, "matplotlib")
    os.makedirs(matplotlib_dir, exist_ok=True)
    os.environ.setdefault("XDG_CACHE_HOME", cache_dir)
    os.environ.setdefault("MPLCONFIGDIR", matplotlib_dir)
    warnings.filterwarnings("ignore", category=SyntaxWarning)

    try:
        from sytogen.scripts import legacy_sytogen
    except ModuleNotFoundError as exc:
        missing = exc.name or "a legacy SyToGen dependency"
        print("IMPORT ERROR:", repr(exc))
        raise SyToGenPipelineError(
            f"Legacy SyToGen dependency is not installed: {missing}"
        ) from exc

    return legacy_sytogen


def run_legacy_estimator(paths, params):
    input_sequence = _first_existing(paths, "genbank")
    input_rm_systems = _first_existing(paths, "motif_table")
    input_strain_genome = _first_existing(paths, "genome", "strain_genome", "genbank")
    missing = []
    print("PATHS RECEIVED:", paths)
    if not input_sequence:
        missing.append("GenBank construct file")
    if not input_rm_systems:
        missing.append("RM motif table")
    if not input_strain_genome:
        missing.append("strain genome GenBank file")
    if missing:
        raise SyToGenPipelineError(
            "Legacy SyToGen requires: " + ", ".join(missing))

    legacy = _legacy_module()
    output_dir = params["output_dir"]

    args = SimpleNamespace(
        input_sequence=input_sequence,
        input_rm_systems=input_rm_systems,
        input_strain_genome=input_strain_genome,
        codon_table=int(params.get("codon_table") or 11),
        output_folder=output_dir,
        verbose=False,
        NCORES=int(params.get("ncores") or 1),
    )

    legacy.run_estimator(args)

    return {
        "input_sequence": input_sequence,
        "input_rm_systems": input_rm_systems,
        "output_dir": output_dir,
    }


def run_legacy_candidate_builder(step1, params):
    legacy = _legacy_module()
    output_dir = step1["output_dir"]

    args = SimpleNamespace(
        input_sequence=step1["input_sequence"],
        input_rm_systems=step1["input_rm_systems"],
        synpl_folder=output_dir,
        output_folder=output_dir,
        auto=int(params.get("auto") or 1),
        verbose=False,
    )

    legacy.candidate_builder(args)

    fasta_path = os.path.join(output_dir, "candidate_syngenic_sequence__auto__.fa")
    genbank_path = os.path.join(output_dir, "candidate_syngenic_sequence__auto__.gbk")

    for path in (fasta_path, genbank_path):
        if os.path.exists(path):
            with open(path) as handle:
                return handle.read()

    raise SyToGenPipelineError("Legacy SyToGen did not produce a candidate sequence")
