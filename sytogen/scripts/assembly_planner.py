"""
assembly_planner.py

Genome-aware Gibson Assembly planning for SyToGen.

Features
--------
* Soft-avoids engineered edit locations, coding regions, and protected
  regions when choosing fragment boundaries (penalized during boundary
  search, not hard-excluded — see warnings below for when that
  compromise actually happened)
* Optimizes overlap placement
* Produces assembly QA metrics, including human-readable warnings when
  boundary/overlap selection had to compromise
* Supports circular plasmids
* Designs forward/reverse PCR primers for each fragment, with the shared
  Gibson homology overlap as the 5' tail and a Tm-targeted gene-specific
  3' annealing region
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional
from .sequence_utils import reverse_complement


# =============================================================================
# Configuration
# =============================================================================

DEFAULT_FRAGMENT_SIZE = 1500
DEFAULT_OVERLAP_LENGTH = 35

BOUNDARY_SEARCH_WINDOW = 250

TARGET_GC = 50.0
TARGET_TM = 65.0

# Primer annealing region (the gene-specific 3' portion, NOT counting the
# Gibson homology tail — the tail doesn't anneal to the template on the
# first PCR cycle, so it's excluded from these Tm/length targets, which is
# standard practice for Gibson/NEBuilder-style primer design).
PRIMER_MIN_ANNEAL = 18
PRIMER_MAX_ANNEAL = 30
PRIMER_TARGET_TM = 60.0

# --- Nearest-neighbor Tm (SantaLucia 1998 unified parameters) --------------
# Replaces the simple Wallace rule (2*(A+T) + 4*(G+C)) that was used
# everywhere Tm mattered — Gibson overlap scoring and primer annealing
# design. Wallace's rule was only ever a rough estimate for very short
# (<14bp) oligos and ignores sequence context, salt, and oligo
# concentration entirely. Nearest-neighbor thermodynamics accounts for
# stacking energy between adjacent base pairs, which is what actually
# governs duplex stability.
#
# ΔH (kcal/mol) / ΔS (cal/(mol·K)) per nearest-neighbor step, indexed by
# the top-strand 5'->3' dinucleotide. Source: SantaLucia, J. PNAS 1998,
# 95(4):1460-1465, Table 1 ("unified" parameters).
NN_THERMO_PARAMS = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
# Helix-initiation penalties, applied once per terminal base (5' and 3' ends).
NN_INIT_AT = (2.3, 4.1)     # ΔH, ΔS for a terminal A or T
NN_INIT_GC = (0.1, -2.8)    # ΔH, ΔS for a terminal G or C

GAS_CONSTANT_CAL = 1.987  # cal/(mol*K)

# Typical PCR/Gibson defaults — 250nM primer/oligo, 50mM monovalent cation
# (close to standard Taq/Gibson Master Mix buffer conditions). Exposed as
# module constants so they're easy to tune for a specific protocol.
DEFAULT_OLIGO_CONC_M = 250e-9
DEFAULT_MONOVALENT_CONC_M = 50e-3


# =============================================================================
# Data Models
# =============================================================================

@dataclass
class Overlap:
    start: int
    end: int
    sequence: str
    tm: float
    gc_percent: float
    score: float


@dataclass
class Primer:
    name: str
    sequence: str
    tm: float


@dataclass
class AssemblyFragment:
    name: str
    start: int
    end: int
    sequence: str

    overlap_left: Optional[Overlap] = None
    overlap_right: Optional[Overlap] = None

    forward_primer: Optional[Primer] = None
    reverse_primer: Optional[Primer] = None


@dataclass
class AssemblyPlan:
    fragments: List[AssemblyFragment]
    warnings: List[str]
    assembly_score: float
    overlap_length: int
    target_fragment_size: int


# =============================================================================
# Utilities
# =============================================================================

_RC = str.maketrans(
    "ACGTacgt",
    "TGCAtgca"
)


def reverse_complement(seq: str) -> str:
    if not isinstance(seq, str):
        raise TypeError("Sequence must be provided as a string.")
    if any(base not in "ACGTacgt" for base in seq):
        raise ValueError("Sequence contains non-canonical bases.")
    return seq.translate(_RC)[::-1]


def gc_percent(seq: str) -> float:
    seq = seq.upper()

    if not seq:
        return 0.0

    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)


def wallace_tm(seq: str) -> float:
    """
    The old, crude estimate: 2*(A+T) + 4*(G+C). Kept only as a fallback
    for sequences nearest_neighbor_tm() can't handle (too short for any
    NN step, or containing a non-ACGT character) — never used as the
    primary Tm model anymore.
    """
    seq = seq.upper()

    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")

    return (2 * at) + (4 * gc)


def nearest_neighbor_tm(
    seq: str,
    oligo_conc_m: float = DEFAULT_OLIGO_CONC_M,
    monovalent_conc_m: float = DEFAULT_MONOVALENT_CONC_M,
) -> float:
    """
    Melting temperature via nearest-neighbor thermodynamics (SantaLucia
    1998 unified parameters) with a salt correction for monovalent cation
    concentration — the standard model tools like Primer3/IDT OligoAnalyzer
    use, and a meaningfully better estimate than the Wallace rule for any
    oligo long enough to matter for PCR or Gibson overlap design (which is
    always the case here: PRIMER_MIN_ANNEAL=18, DEFAULT_OVERLAP_LENGTH=35).

    Falls back to wallace_tm() for a sequence with no valid nearest-neighbor
    step (shorter than 2bp, or containing an ambiguous/non-ACGT base) —
    that's always an edge case here, never the normal path.
    """
    seq = seq.upper()

    if len(seq) < 2 or any(base not in "ACGT" for base in seq):
        return wallace_tm(seq)

    delta_h = 0.0
    delta_s = 0.0

    for end_base in (seq[0], seq[-1]):
        h, s = NN_INIT_AT if end_base in "AT" else NN_INIT_GC
        delta_h += h
        delta_s += s

    for i in range(len(seq) - 1):
        step = seq[i:i + 2]
        if step not in NN_THERMO_PARAMS:
            return wallace_tm(seq)  # shouldn't happen once we've checked ACGT-only, but stay safe
        h, s = NN_THERMO_PARAMS[step]
        delta_h += h
        delta_s += s

    # Tm in Kelvin for a non-self-complementary duplex (the standard
    # C_T/4 convention — two different strands, one at each end of a
    # linear PCR product / Gibson overlap, never present at equal molar
    # self-complementary concentration).
    tm_kelvin = (delta_h * 1000.0) / (delta_s + GAS_CONSTANT_CAL * math.log(oligo_conc_m / 4))

    # Salt correction (Owczarzy et al. 2004 / SantaLucia's standard
    # 1/Tm-form correction) for monovalent cation concentration —
    # the raw NN parameters above are calibrated for 1M NaCl, far above
    # any real PCR/Gibson buffer.
    gc_fraction = (seq.count("G") + seq.count("C")) / len(seq)
    ln_na = math.log(monovalent_conc_m)
    inv_tm_salt_corrected = (
        (1.0 / tm_kelvin)
        + (4.29 * gc_fraction - 3.95) * 1e-5 * ln_na
        + 9.4e-6 * (ln_na ** 2)
    )
    tm_kelvin_corrected = 1.0 / inv_tm_salt_corrected

    return tm_kelvin_corrected - 273.15


def longest_homopolymer(seq: str) -> int:

    if not seq:
        return 0

    longest = 1
    current = 1

    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current += 1
        else:
            current = 1

        longest = max(longest, current)

    return longest


# =============================================================================
# Decision Matrix Integration
# =============================================================================

def collect_edit_positions(
    decision_matrix
) -> set:

    edited = set()

    for row in decision_matrix:

        if not row.get("chosen"):
            continue

        pos = row.get("edit_position")

        if pos in ("", None):
            continue

        edited.add(int(pos))

    return edited


# =============================================================================
# Genome Helpers
# =============================================================================

def is_in_gene(
    genome,
    position: int
) -> bool:

    gene = genome.find_gene(position)

    return gene is not None


def is_protected(
    genome,
    position: int
) -> bool:

    for region in getattr(genome, "protected_regions", []):

        if region.start <= position <= region.end:
            return True

    return False


def nearest_edit_distance(
    position: int,
    edits: set[int]
) -> int:

    if not edits:
        return 999999

    return min(abs(position - x) for x in edits)


# =============================================================================
# Overlap Scoring
# =============================================================================

def overlap_score(
    sequence: str
) -> tuple[float, float, float]:

    gc = gc_percent(sequence)
    tm = nearest_neighbor_tm(sequence)

    score = 100.0

    score -= abs(gc - TARGET_GC)

    score -= abs(tm - TARGET_TM)

    if longest_homopolymer(sequence) > 5:
        score -= 50

    if sequence.count(sequence[:8]) > 1:
        score -= 15

    return score, tm, gc


# =============================================================================
# Boundary Selection
# =============================================================================

def boundary_score(
    sequence: str,
    position: int,
    genome,
    edited_positions: set[int],
    overlap_length: int,
) -> float:

    score = 0.0

    distance = nearest_edit_distance(
        position,
        edited_positions
    )

    score += min(distance, 200)

    if is_in_gene(genome, position):
        score -= 1000

    if is_protected(genome, position):
        score -= 1000

    overlap_start = max(
        0,
        position - overlap_length // 2
    )

    overlap_end = min(
        len(sequence),
        overlap_start + overlap_length
    )

    overlap_seq = sequence[
        overlap_start:overlap_end
    ]

    overlap_quality, _, _ = overlap_score(
        overlap_seq
    )

    score += overlap_quality

    return score


def choose_boundary(
    sequence: str,
    approximate_boundary: int,
    genome,
    edited_positions: set[int],
    overlap_length: int,
) -> int:

    best_position = approximate_boundary
    best_score = float("-inf")

    start = max(
        overlap_length,
        approximate_boundary - BOUNDARY_SEARCH_WINDOW
    )

    end = min(
        len(sequence) - overlap_length,
        approximate_boundary + BOUNDARY_SEARCH_WINDOW
    )

    for pos in range(start, end):

        score = boundary_score(
            sequence,
            pos,
            genome,
            edited_positions,
            overlap_length,
        )

        if score > best_score:

            best_score = score
            best_position = pos

    return best_position


# =============================================================================
# Fragment Generation
# =============================================================================

def create_overlap(
    sequence: str,
    boundary: int,
    overlap_length: int
) -> Overlap:

    start = boundary - (overlap_length // 2)

    end = start + overlap_length

    overlap_seq = sequence[start:end]

    score, tm, gc = overlap_score(
        overlap_seq
    )

    return Overlap(
        start=start,
        end=end,
        sequence=overlap_seq,
        tm=tm,
        gc_percent=gc,
        score=score,
    )


def create_overlap_at_origin(
    sequence: str,
    overlap_length: int
) -> Overlap:
    """
    Same as create_overlap(), but wrap-aware — for the boundary=0 junction
    on a circular molecule, where the "boundary" isn't really an edge at
    all, just an arbitrary reference point in a continuous loop. Half the
    overlap window sits at the end of `sequence` and half at the start.

    Overlap.start / Overlap.end are left un-clamped here (start can be
    negative, end can exceed len(sequence)) the same way Motif.start/end
    are for an origin-spanning motif — every other place in the codebase
    that positions primers off an Overlap already normalizes coordinates
    the same way genome_model.GenomeModel's topology engine does.
    """
    length = len(sequence)
    start = -(overlap_length // 2)
    end = start + overlap_length

    overlap_seq = sequence[start:] + sequence[:end]

    score, tm, gc = overlap_score(
        overlap_seq
    )

    return Overlap(
        start=start,
        end=end,
        sequence=overlap_seq,
        tm=tm,
        gc_percent=gc,
        score=score,
    )


def _extract_window(sequence: str, start: int, length: int, circular: bool) -> str:
    """
    `length` bases starting at `start`, extending forward (5'->3' on the
    plus strand). Wraps for circular topology; truncates at the physical
    end of the sequence for linear topology (start is normalized via
    modulo first, same convention used for Motif/Overlap positions
    elsewhere in this codebase — but only when it's actually out of
    range, since len(sequence) % len(sequence) == 0 would otherwise wrap
    a valid boundary position back to the start).
    """
    L = len(sequence)
    if not (0 <= start < L):
        start %= L
    end = start + length

    if end <= L:
        return sequence[start:end]
    if circular:
        return sequence[start:] + sequence[:end - L]
    return sequence[start:L]


def _extract_window_reverse(sequence: str, end_exclusive: int, length: int, circular: bool) -> str:
    """
    The `length` bases immediately preceding `end_exclusive` (still on the
    plus strand, 5'->3' — the caller reverse-complements afterward to get
    an actual reverse-primer annealing region). Wraps for circular
    topology; truncates at the physical start of the sequence for linear.

    end_exclusive is only normalized via modulo when it's actually out of
    the valid [0, L] range — len(sequence) is itself a legitimate
    exclusive upper bound (the true end of a linear fragment) and must
    NOT be wrapped to 0, or the last fragment's reverse primer would be
    built from an empty annealing region.
    """
    L = len(sequence)
    if not (0 <= end_exclusive <= L):
        end_exclusive %= L
    start = end_exclusive - length

    if start >= 0:
        return sequence[start:end_exclusive]
    if circular:
        return sequence[start % L:] + sequence[:end_exclusive]
    return sequence[0:end_exclusive]


def _grow_anneal_region(extract_fn, min_len=PRIMER_MIN_ANNEAL, max_len=PRIMER_MAX_ANNEAL, target_tm=PRIMER_TARGET_TM):
    """
    Grow an annealing region one base at a time from `min_len` until its
    nearest-neighbor Tm reaches `target_tm` or `max_len` is hit —
    whichever comes first. `extract_fn(n)` returns the n-base candidate
    region.
    """
    length = min_len
    region = extract_fn(length)
    tm = nearest_neighbor_tm(region)

    while tm < target_tm and length < max_len:
        length += 1
        region = extract_fn(length)
        tm = nearest_neighbor_tm(region)

    return region, tm


def design_primers_for_fragment(sequence: str, fragment: AssemblyFragment, circular: bool):
    """
    Build the forward/reverse PCR primer pair for one fragment.

    Each primer = [Gibson homology tail, if this side borders another
    fragment] + [gene-specific annealing region, grown to PRIMER_TARGET_TM].
    The tail is the exact overlap sequence shared with the neighboring
    fragment (so both fragments' amplicons carry the identical homology
    arm), and the annealing region starts right where that tail ends, so
    the two never overlap each other.

    At a true linear terminus (no overlap on that side), the tail is
    empty and the annealing region just starts at the fragment's own edge.
    """
    L = len(sequence)

    # ---- forward primer ----
    if fragment.overlap_left is not None:
        fwd_tail = fragment.overlap_left.sequence
        anneal_start = fragment.overlap_left.end
    else:
        fwd_tail = ""
        anneal_start = fragment.start

    fwd_anneal, fwd_tm = _grow_anneal_region(
        lambda n: _extract_window(sequence, anneal_start, n, circular)
    )

    forward_primer = Primer(
        name=f"{fragment.name}_F",
        sequence=fwd_tail + fwd_anneal,
        tm=fwd_tm,
    )

    # ---- reverse primer ----
    if fragment.overlap_right is not None:
        rev_tail_plus_strand = fragment.overlap_right.sequence
        anneal_end = fragment.overlap_right.start
    else:
        rev_tail_plus_strand = ""
        anneal_end = fragment.end

    rev_anneal_plus_strand, rev_tm = _grow_anneal_region(
        lambda n: _extract_window_reverse(sequence, anneal_end, n, circular)
    )

    reverse_primer = Primer(
        name=f"{fragment.name}_R",
        sequence=reverse_complement(rev_tail_plus_strand) + reverse_complement(rev_anneal_plus_strand),
        tm=rev_tm,
    )

    return forward_primer, reverse_primer


def generate_assembly_warnings(fragments, genome, edited_positions, overlap_length, assembly_score):
    """
    Surface the cases where boundary/overlap selection had to compromise.

    boundary_score() only *penalizes* landing in a gene, a protected
    region, or near an edited position — it never hard-excludes those
    outcomes, so in a gene-dense region a boundary can still end up
    somewhere non-ideal if every candidate position in the search window
    was bad. This function doesn't change that behavior; it just makes
    sure it's visible in the output instead of shipping silently.

    Only checks fragment sides that actually have an overlap (i.e. were
    chosen by choose_boundary/create_overlap) — a linear construct's true
    termini aren't a "boundary" that could have been placed better, so
    they're skipped.
    """
    warnings = []

    for frag in fragments:
        for side, pos, overlap in (
            ("start", frag.start, frag.overlap_left),
            ("end",   frag.end,   frag.overlap_right),
        ):
            if overlap is None:
                continue

            if is_in_gene(genome, pos):
                warnings.append(
                    f"{frag.name}: {side} boundary (position {pos}) falls inside a gene — "
                    f"no clean intergenic cut site was found within the search window."
                )
            elif is_protected(genome, pos):
                warnings.append(
                    f"{frag.name}: {side} boundary (position {pos}) falls inside a "
                    f"protected region — no clean cut site was found within the search window."
                )

            if edited_positions:
                distance = nearest_edit_distance(pos, edited_positions)
                if distance < overlap_length:
                    warnings.append(
                        f"{frag.name}: {side} boundary (position {pos}) is only {distance}bp "
                        f"from a position SyToGen edited — Sanger verification of that edit "
                        f"may be less reliable this close to a primer/junction site."
                    )

        for side, overlap in (("left", frag.overlap_left), ("right", frag.overlap_right)):
            if overlap is not None and overlap.score < 0:
                warnings.append(
                    f"{frag.name}: {side} overlap (Tm={overlap.tm:.1f}C, "
                    f"GC={overlap.gc_percent:.1f}%, score={overlap.score:.1f}) is a poor "
                    f"match for the target Tm/GC — consider a different overlap_length "
                    f"for this construct."
                )

        if len(frag.sequence) < overlap_length * 2:
            warnings.append(
                f"{frag.name}: fragment is only {len(frag.sequence)}bp — close to or below "
                f"twice the overlap length ({overlap_length}bp), so its overlaps may span "
                f"most of the fragment."
            )

    if assembly_score < 0:
        warnings.append(
            f"Overall assembly score is low ({assembly_score:.1f}) — one or more overlaps "
            f"had to compromise on Tm/GC/homopolymer quality; see per-fragment warnings above."
        )

    return warnings


def fragment_sequence(
    sequence: str,
    genome,
    decision_matrix,
    fragment_size=DEFAULT_FRAGMENT_SIZE,
    overlap_length=DEFAULT_OVERLAP_LENGTH,
    topology="circular",
):

    edits = collect_edit_positions(
        decision_matrix
    )

    boundaries = []

    current = fragment_size

    while current < len(sequence):

        boundary = choose_boundary(
            sequence,
            current,
            genome,
            edits,
            overlap_length,
        )

        boundaries.append(boundary)

        current += fragment_size

    # ensure start and end boundaries
    if not boundaries or boundaries[0] != 0:
        boundaries.insert(0, 0)

    if boundaries[-1] != len(sequence):
        boundaries.append(len(sequence))

    fragments: List[AssemblyFragment] = []

    for i in range(len(boundaries) - 1):
        start = boundaries[i]
        end = boundaries[i + 1]
        seq = sequence[start:end]

        frag = AssemblyFragment(
            name=f"frag{i+1}",
            start=start,
            end=end,
            sequence=seq,
        )

        # overlaps
        if i > 0:
            frag.overlap_left = create_overlap(sequence, start, overlap_length)

        if i < len(boundaries) - 2:
            frag.overlap_right = create_overlap(sequence, end, overlap_length)

        fragments.append(frag)

    # For a circular molecule, position 0 isn't a real edge — the last
    # fragment needs to overlap back into the first one to close the ring.
    # Fragment boundaries never place a cut exactly at the origin unless
    # there's only one fragment total, so this is safe to add unconditionally.
    if topology == "circular" and len(fragments) > 1:
        origin_overlap = create_overlap_at_origin(sequence, overlap_length)
        fragments[0].overlap_left = origin_overlap
        fragments[-1].overlap_right = origin_overlap

    # simple assembly score: average overlap score
    scores = [
        o.score
        for f in fragments
        for o in (f.overlap_left, f.overlap_right)
        if o is not None
    ]

    assembly_score = sum(scores) / len(scores) if scores else 0.0

    # Primer design — every fragment's forward/reverse primer pair, with
    # the shared homology overlap as the 5' tail (empty at the true
    # termini of a linear construct, where there's nothing to join to).
    for frag in fragments:
        frag.forward_primer, frag.reverse_primer = design_primers_for_fragment(
            sequence, frag, topology == "circular"
        )

    warnings = generate_assembly_warnings(
        fragments, genome, edits, overlap_length, assembly_score
    )

    return AssemblyPlan(
        fragments=fragments,
        warnings=warnings,
        assembly_score=assembly_score,
        overlap_length=overlap_length,
        target_fragment_size=fragment_size,
    )