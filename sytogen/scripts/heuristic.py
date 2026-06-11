# Heuristic optimization of a sequence based on motif positions and codon synonyms.
# This is a simple greedy approach that iteratively replaces codons in motif positions 
# with their synonyms, choosing the one that reduces the number of remaining motifs the most.


import copy
import math
import re
from Bio.Seq import reverse_complement
from Bio.Seq import translate
from sytogen.scripts.pipeline import count_motifs
from sytogen.scripts.pipeline import run_legacy_candidate_builder

# This is a simple heuristic optimization that iteratively replaces codons in 
# motif positions with their synonyms, choosing the one that reduces the number 
# of remaining motifs the most. 

def score_sequence(seq, cds_repr, motifs, weights):
    """
    Global scoring function.
    Lower = better.
    """
    # motif penalty = number of motifs in the sequence
    seq_str = "".join(seq)

    motif_penalty = count_motifs(seq_str)

    codon_penalty = 0
    edit_penalty = 0

    # codon penalty = sum of -log(usage) for each codon used
    for pos, info in cds_repr.items():
        if info["type"] != "codon":
            continue

        original = str(info["original"])
        current = seq[pos:pos+3]

        if "".join(current) != original:
            edit_penalty += 1

        for codon, usage in info["synonims"].items():
            if "".join(current) == codon:
                codon_penalty += -math.log(usage + 1e-6)
                break
    # combine penalties with weights
    return (
        weights["motif"] * motif_penalty +
        weights["codon"] * codon_penalty +
        weights["edit"] * edit_penalty
    )

# This is a simple heuristic optimization that iteratively replaces codons in 
# motif positions with their synonyms, choosing the one that reduces the number 
# of remaining motifs the most.

def heuristic_optimize(sequence, motifs, cds_repr, beam_width=5):

    weights = {
        "motif": 10.0,
        "codon": 1.0,
        "edit": 0.5
    }

    initial = list(sequence)

    # each state = (sequence, score)
    beam = [(initial, score_sequence(initial, cds_repr, motifs, weights))]

    for motif in motifs:
        positions = motif["positions"]

        candidates = []

        for seq, _ in beam:

            for pos in positions:
                info = cds_repr.get(pos)
                if not info or info["type"] != "codon":
                    continue

                original = str(info["original"])

                for codon in info["synonims"]:

                    if codon == original:
                        continue

                    new_seq = seq.copy()
                    new_seq[pos:pos+3] = list(codon)

                    score = score_sequence(new_seq, cds_repr, motifs, weights)

                    candidates.append((new_seq, score))

        # keep best k candidates
        candidates.sort(key=lambda x: x[1])
        beam = candidates[:beam_width]

        if not beam:
            break

    # return best sequence
    best_seq = min(beam, key=lambda x: x[1])[0]
    return "".join(best_seq)

def check_synonymous(new_seq, original_seq, cds_repr):
    # loop through CDS regions
    for start, end, strand in cds_repr:
        orig = original_seq[start:end]
        new = new_seq[start:end]

        if strand == "-":
            orig = reverse_complement(orig)
            new = reverse_complement(new)

        if translate(orig) != translate(new):
            return False

    return True

def compute_codon_bias(sequence, cds_repr, codon_usage):
    penalty = 0

    for start, end, strand in cds_repr:
        seq = sequence[start:end]

        if strand == "-":
            seq = reverse_complement(seq)

        codons = re.findall(r'.{3}', seq)

        for codon in codons:
            if codon in codon_usage:
                penalty += codon_usage[codon]["Inverse_ranking"]

    return penalty

def run_step2(sequence, 
              original_sequence, 
              motifs, 
              cds_repr, 
              codon_usage, 
              weights=None):
    
    weights = weights or {
        "motif": 1.0,
        "codon": 0.1,
        "modifications": 0.01,
        "nonsynonymous_penalty": 1000}

    # --- 1. Motif count ---
    motif_count = count_motifs(sequence, motifs)

    # --- 2. Synonymous check ---
    is_synonymous = check_synonymous(sequence, original_sequence, cds_repr)

    nonsyn_penalty = 0 if is_synonymous else weights["nonsynonymous_penalty"]

    # --- 3. Codon bias ---
    codon_penalty = compute_codon_bias(sequence, cds_repr, codon_usage)

    # --- 4. Number of changes ---
    modification_count = sum(
        1 for a, b in zip(original_sequence, sequence) if a != b)

    # --- 5. Final score (lower is better) ---
    score = (weights["motif"] * motif_count
        + weights["codon"] * codon_penalty
        + weights["modifications"] * modification_count
        + nonsyn_penalty)

    return {"score": score,
        "motif_count": motif_count,
        "codon_penalty": codon_penalty,
        "modifications": modification_count,
        "is_synonymous": is_synonymous}