import pulp
import math
from sytogen.scripts.pipeline import run_legacy_estimator

def solve_ilp_for_segment(seq, cds_repr, weights, forced=None, excluded=None):
    model = pulp.LpProblem("sytogen_ilp", pulp.LpMinimize)

    x = {}

    # build variables
    for pos, info in cds_repr.items():
        if info["type"] != "codon":
            continue

        x[pos] = {}

        for codon in info["synonims"].keys():
            name = f"x_{pos}_{codon}"
            x[pos][codon] = pulp.LpVariable(name, cat="Binary")

    # 1. one codon constraint
    for pos in x:
        model += pulp.lpSum(x[pos][c] for c in x[pos]) == 1

    # 2. forced edits
    if forced:
        for pos, codon in forced.items():
            if pos in x and codon in x[pos]:
                model += x[pos][codon] == 1

    # 3. excluded edits
    if excluded:
        for pos, codons in excluded.items():
            for codon in codons:
                if pos in x and codon in x[pos]:
                    model += x[pos][codon] == 0

    # objective
    obj = []

    for pos in x:
        info = cds_repr[pos]

        for codon, var in x[pos].items():

            usage = info["synonims"][codon]

            codon_cost = -math.log(usage + 1e-6)

            obj.append(weights["codon"] * codon_cost * var)

    model += pulp.lpSum(obj)

    model.solve(pulp.PULP_CBC_CMD(msg=False))

    # extract solution
    result = {}

    for pos in x:
        for codon, var in x[pos].items():
            if pulp.value(var) == 1:
                result[pos] = codon

    return result

import copy


def score_full_sequence(seq, motifs):
    return count_motifs("".join(seq))


def apply_ilp_solution(seq, solution):
    new_seq = seq.copy()

    for pos, codon in solution.items():
        new_seq[pos:pos+3] = list(codon)

    return new_seq


def beam_ilp_optimize(sequence, motifs, cds_repr,
                      beam_width=5,
                      ilp_chunk_size=10):

    weights = {
        "motif": 10.0,
        "codon": 1.0,
        "edit": 0.5
    }

    initial = list(sequence)
    beam = [(initial, score_full_sequence(initial, motifs))]

    codon_positions = [p for p, v in cds_repr.items() if v["type"] == "codon"]

    for start in range(0, len(codon_positions), ilp_chunk_size):

        chunk = codon_positions[start:start+ilp_chunk_size]

        new_candidates = []

        for seq, _ in beam:

            # build cds subset
            sub_repr = {p: cds_repr[p] for p in chunk}

            solution = solve_ilp_for_segment(
                seq,
                sub_repr,
                weights,
                forced=None,
                excluded=None
            )

            new_seq = apply_ilp_solution(seq, solution)

            score = score_full_sequence(new_seq, motifs)

            new_candidates.append((new_seq, score))

        # beam prune
        new_candidates.sort(key=lambda x: x[1])
        beam = new_candidates[:beam_width]

    return min(beam, key=lambda x: x[1])[0]


def run_step1(paths, params):
    # print(locals())
    return paths
