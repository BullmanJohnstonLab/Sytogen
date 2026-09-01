"""Constraint-based DNA sequence optimization utilities."""

import re

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

from .dna_analysis import find_dispersed_repeats, find_tandem_repeats
from .pattern_library import PatternLibrary, get_default_library, iupac_to_regex
from .sequence_utils import gc_percent, longest_homopolymer


@dataclass
class OptimizationConstraints:
    """Constraints applied by :func:`optimize_sequence`."""

    forbidden_patterns: List[str] = field(default_factory=list)
    min_gc: Optional[float] = None
    max_gc: Optional[float] = None
    max_homopolymer: Optional[int] = None
    avoid_tandem_repeats: bool = True
    avoid_dispersed_repeats: bool = True
    tandem_repeat_min_length: int = 4
    dispersed_repeat_min_length: int = 8

    @classmethod
    def from_dict(cls, values: Optional[Dict] = None) -> "OptimizationConstraints":
        values = values or {}
        allowed = {
            "forbidden_patterns", "min_gc", "max_gc", "max_homopolymer",
            "avoid_tandem_repeats", "avoid_dispersed_repeats",
            "tandem_repeat_min_length", "dispersed_repeat_min_length",
        }
        unknown = set(values) - allowed
        if unknown:
            raise ValueError(f"Unknown constraint(s): {', '.join(sorted(unknown))}")
        constraints = cls(**values)
        if constraints.min_gc is not None and constraints.max_gc is not None:
            if constraints.min_gc > constraints.max_gc:
                raise ValueError("min_gc cannot exceed max_gc")
        if constraints.max_homopolymer is not None and constraints.max_homopolymer < 1:
            raise ValueError("max_homopolymer must be at least 1")
        return constraints


def _forbidden_matches(sequence: str, patterns: Iterable[str], library: PatternLibrary) -> List[Tuple[str, int, int]]:
    matches = []
    for pattern in patterns:
        spec = library.get_pattern(pattern)
        regex = spec.compile() if spec else re.compile(iupac_to_regex(pattern), re.IGNORECASE)
        matches.extend((pattern, hit.start(), hit.end()) for hit in regex.finditer(sequence))
    return matches


def constraint_violations(
    sequence: str,
    constraints: OptimizationConstraints,
    library: Optional[PatternLibrary] = None,
) -> Dict[str, int]:
    """Return an auditable violation count for each configured constraint."""
    sequence = sequence.upper()
    library = library or get_default_library()
    violations = {
        "forbidden_patterns": len(_forbidden_matches(sequence, constraints.forbidden_patterns, library)),
        "gc_bounds": 0,
        "homopolymers": 0,
        "tandem_repeats": 0,
        "dispersed_repeats": 0,
    }
    gc = gc_percent(sequence)
    if constraints.min_gc is not None and gc < constraints.min_gc:
        violations["gc_bounds"] = 1
    if constraints.max_gc is not None and gc > constraints.max_gc:
        violations["gc_bounds"] = 1
    if constraints.max_homopolymer is not None:
        violations["homopolymers"] = int(longest_homopolymer(sequence) > constraints.max_homopolymer)
    if constraints.avoid_tandem_repeats:
        violations["tandem_repeats"] = len(find_tandem_repeats(sequence, constraints.tandem_repeat_min_length))
    if constraints.avoid_dispersed_repeats:
        violations["dispersed_repeats"] = len(find_dispersed_repeats(sequence, constraints.dispersed_repeat_min_length))
    return violations


def optimize_sequence(
    sequence: str,
    constraints: OptimizationConstraints,
    max_edits: int = 100,
    library: Optional[PatternLibrary] = None,
) -> Dict:
    """Greedily apply substitutions that reduce the total constraint burden.

    This is a general sequence optimizer and does not preserve translated
    protein sequence. Coding-region optimization remains SyToGen's job,
    where synonymous-codon constraints are available.
    """
    sequence = sequence.upper()
    invalid = set(sequence) - set("ACGT")
    if invalid:
        raise ValueError(f"Sequence contains non-canonical bases: {sorted(invalid)}")
    if max_edits < 0:
        raise ValueError("max_edits must be non-negative")

    library = library or get_default_library()
    optimized = sequence
    edits = []
    before = constraint_violations(optimized, constraints, library)

    for _ in range(max_edits):
        current = constraint_violations(optimized, constraints, library)
        current_total = sum(current.values())
        if current_total == 0:
            break

        best = None
        for position, old_base in enumerate(optimized):
            for new_base in "ACGT":
                if new_base == old_base:
                    continue
                candidate = optimized[:position] + new_base + optimized[position + 1:]
                candidate_violations = constraint_violations(candidate, constraints, library)
                candidate_total = sum(candidate_violations.values())
                if candidate_total < current_total:
                    choice = (candidate_total, position, old_base, new_base, candidate, candidate_violations)
                    if best is None or choice[:3] < best[:3]:
                        best = choice
        if best is None:
            break
        _, position, old_base, new_base, optimized, _ = best
        edits.append({"position": position, "original": old_base, "edited": new_base})

    after = constraint_violations(optimized, constraints, library)
    return {
        "original_sequence": sequence,
        "optimized_sequence": optimized,
        "edits": edits,
        "edits_applied": len(edits),
        "violations_before": before,
        "violations_after": after,
        "resolved": sum(after.values()) == 0,
    }
