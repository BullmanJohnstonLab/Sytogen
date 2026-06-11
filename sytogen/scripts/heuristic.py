from sytogen.scripts.genome_model import GenomeModel, Mutation, evaluate_mutation


# =========================================================
# Heuristic optimizer
# =========================================================
def heuristic_optimize(sequence, genome_model, cds_repr, beam_width=5):

    beam = [{
        "sequence": sequence,
        "score": 0,
        "mutations": []
    }]

    for pos, info in cds_repr.items():

        if info.get("type") != "codon":
            continue

        candidates = []

        for candidate in beam:
            seq = candidate["sequence"]
            original_codon = seq[pos:pos+3]
            synonyms = info.get("synonyms", [])

            for codon in synonyms:

                if codon == original_codon:
                    continue

                mutation = Mutation(
                    position=pos,
                    old=original_codon,
                    new=codon
                )

                # Build per-candidate model
                temp_model = GenomeModel(
                    sequence=seq,
                    topology=genome_model.topology,
                    genes=genome_model.genes,
                    motifs=genome_model.motifs,
                    protected_regions=genome_model.protected_regions,
                    codon_usage=genome_model.codon_usage
                )

                result = evaluate_mutation(temp_model, mutation)

                if not result["valid"]:
                    continue

                # Apply mutation
                new_seq = list(seq)
                new_seq[pos:pos+3] = list(codon)
                new_seq = "".join(new_seq)

                score = (
                    candidate["score"]
                    - 10 * result["destroyed"]
                    + result["edits"]
                )

                candidates.append({
                    "sequence": new_seq,
                    "score": score,
                    "mutations": candidate["mutations"] + [mutation]
                })

        # beam pruning
        if candidates:
            candidates.sort(key=lambda x: x["score"])
            beam = candidates[:beam_width]

    return beam[0]["sequence"]


# =========================================================
# Pipeline entry point (API uses this)
# =========================================================
def run_step2(
    sequence,
    original_sequence,
    motifs,
    cds_repr,
    codon_usage,
    weights=None,
):
    """
    Replacement for legacy Step 2.
    """

    print("=== NEW HEURISTIC ENGINE RUNNING ===")
    # --- CONVERT MOTIFS TO OBJECTS ---
    converted_motifs = []

    for m in motifs:
        try:
            converted_motifs.append(
                Motif(
                    motif=m["motif"],
                    start=m["start"],
                    end=m["end"],
                    strand=m.get("strand", "+")
                )
            )
        except Exception as e:
            print("Motif conversion error:", m, e)

    # --- BUILD GENOME MODEL ---
    genome_model = GenomeModel(
        sequence=sequence,
        topology="circular",  # TODO: make configurable
        genes=[],
        motifs=converted_motifs,
        protected_regions=[],
        codon_usage=codon_usage,
    )

    # --- RUN OPTIMIZER ---
    optimized_seq = heuristic_optimize(
        sequence=sequence,
        genome_model=genome_model,
        cds_repr=cds_repr,
        beam_width=5,
    )

    return {
        "optimized_sequence": optimized_seq,
        "original_sequence": original_sequence,
        "motifs_remaining": None,
    }