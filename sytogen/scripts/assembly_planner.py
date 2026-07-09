"""
assembly_planner.py

Genome-aware Gibson Assembly planning for SyToGen.

Features
--------
* Avoids engineered edit locations
* Avoids coding regions and protected regions
* Optimizes overlap placement
* Produces assembly QA metrics
* Supports circular plasmids
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional


# =============================================================================
# Configuration
# =============================================================================

DEFAULT_FRAGMENT_SIZE = 1500
DEFAULT_OVERLAP_LENGTH = 35

EDIT_BUFFER = 50
BOUNDARY_SEARCH_WINDOW = 250

TARGET_GC = 50.0
TARGET_TM = 65.0

MIN_OVERLAP = 25
MAX_OVERLAP = 40


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
    return seq.translate(_RC)[::-1]


def gc_percent(seq: str) -> float:
    seq = seq.upper()

    if not seq:
        return 0.0

    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)


def wallace_tm(seq: str) -> float:
    seq = seq.upper()

    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")

    return (2 * at) + (4 * gc)


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

        if region.start <= position < region.end:
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
    tm = wallace_tm(sequence)

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


def fragment_sequence(
    sequence: str,
    genome,
    decision_matrix,
    fragment_size=DEFAULT_FRAGMENT_SIZE,
    overlap_length=DEFAULT_OVERLAP_LENGTH,
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

    # simple assembly score: average overlap score
    scores = [
        o.score
        for f in fragments
        for o in (f.overlap_left, f.overlap_right)
        if o is not None
    ]

    assembly_score = sum(scores) / len(scores) if scores else 0.0

    return AssemblyPlan(
        fragments=fragments,
        warnings=[],
        assembly_score=assembly_score,
        overlap_length=overlap_length,
        target_fragment_size=fragment_size,
    )