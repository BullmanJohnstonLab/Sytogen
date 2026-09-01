"""
DNA sequence analysis utilities for synthetic biology.

Features:
- Repeat detection (tandem and dispersed)
- Secondary structure prediction
- Homopolymer detection
- Sequence quality scoring
"""

import re
from typing import List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class Repeat:
    """Represents a sequence repeat."""
    sequence: str
    positions: List[int]
    repeat_length: int
    copy_count: int
    spacing: Optional[int] = None  # For tandem repeats
    
    def __repr__(self):
        if self.spacing == 0:
            return (f"TandemRepeat({self.sequence!r}, {self.copy_count}x "
                   f"at positions {self.positions}, {self.repeat_length}bp)")
        else:
            return (f"DispersedRepeat({self.sequence!r}, {self.copy_count}x "
                   f"at positions {self.positions})")


@dataclass
class SecondaryStructure:
    """Represents a predicted secondary structure region."""
    start: int
    end: int
    structure_type: str  # "hairpin", "loop", "bulge", etc.
    stability_kcal: float
    sequence: str
    
    def __repr__(self):
        return (f"{self.structure_type} at {self.start}-{self.end} "
               f"(ΔG={self.stability_kcal:.1f} kcal/mol)")


def find_tandem_repeats(sequence: str, min_repeat_length: int = 4, 
                        max_repeat_copies: int = 10) -> List[Repeat]:
    """
    Find tandem repeats (direct repeats at same location).
    
    Args:
        sequence: DNA sequence
        min_repeat_length: Minimum motif length to detect
        max_repeat_copies: Maximum number of copies to allow
        
    Returns:
        List of Repeat objects representing tandem repeats
    """
    sequence = sequence.upper()
    repeats = []
    
    for repeat_len in range(min_repeat_length, len(sequence) // 2):
        for start in range(len(sequence) - repeat_len):
            motif = sequence[start:start + repeat_len]
            
            # Count consecutive copies
            pos = start
            copy_count = 0
            positions = []
            
            while pos + repeat_len <= len(sequence):
                if sequence[pos:pos + repeat_len] == motif:
                    positions.append(pos)
                    copy_count += 1
                    pos += repeat_len
                else:
                    break
            
            # Report if we found multiple copies
            if copy_count >= 2 and copy_count <= max_repeat_copies:
                repeats.append(Repeat(
                    sequence=motif,
                    positions=positions,
                    repeat_length=repeat_len,
                    copy_count=copy_count,
                    spacing=0
                ))
    
    # Remove duplicates (keep longest)
    unique = {}
    for rep in repeats:
        key = (rep.sequence, tuple(rep.positions))
        if key not in unique or rep.repeat_length > unique[key].repeat_length:
            unique[key] = rep
    
    return list(unique.values())


def find_dispersed_repeats(sequence: str, min_repeat_length: int = 8,
                          max_distance: int = 500) -> List[Repeat]:
    """
    Find dispersed repeats (same sequence at different locations).
    Useful for detecting recombination hotspots.
    
    Args:
        sequence: DNA sequence
        min_repeat_length: Minimum motif length
        max_distance: Maximum distance between repeats
        
    Returns:
        List of Repeat objects
    """
    sequence = sequence.upper()
    repeats = []
    
    for repeat_len in range(min_repeat_length, min(len(sequence) // 2, 20)):
        seen = {}
        
        for start in range(len(sequence) - repeat_len + 1):
            motif = sequence[start:start + repeat_len]
            
            # Skip low-complexity (all same base)
            if len(set(motif)) == 1:
                continue
            
            if motif not in seen:
                seen[motif] = []
            seen[motif].append(start)
        
        # Report motifs found at 2+ locations
        for motif, positions in seen.items():
            if len(positions) >= 2:
                # Check if they're tandem (handled by find_tandem_repeats)
                is_tandem = all(
                    positions[i+1] - positions[i] == repeat_len 
                    for i in range(len(positions) - 1)
                )
                if not is_tandem:
                    repeats.append(Repeat(
                        sequence=motif,
                        positions=positions,
                        repeat_length=repeat_len,
                        copy_count=len(positions),
                        spacing=positions[1] - positions[0] if len(positions) > 1 else None
                    ))
    
    return repeats


def simple_secondary_structure_score(sequence: str, window_size: int = 40) -> List[Tuple[int, float]]:
    """
    Simple secondary structure risk scoring based on:
    - High GC content (more stable base pairing)
    - Homopolymer runs (potential for hairpins)
    - Inverted repeats (hairpin formation)
    
    Uses heuristic scoring, not full thermodynamic prediction.
    
    Args:
        sequence: DNA sequence
        window_size: Window for local analysis
        
    Returns:
        List of (position, risk_score) tuples. Higher score = more stable structure
    """
    from .sequence_utils import gc_percent, longest_homopolymer
    
    sequence = sequence.upper()
    scores = []
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        
        # GC content drives stability
        gc = gc_percent(window)
        gc_score = abs(gc - 50) / 50  # Peak at 50% GC
        
        # Homopolymer runs indicate potential structure
        hpoly = longest_homopolymer(window)
        hpoly_score = min(hpoly / 10, 1.0)  # Saturate at 10bp runs
        
        # Check for inverted repeats (hairpin potential)
        inverted_score = 0.0
        half = window_size // 2
        left = window[:half]
        right = window[half:]
        
        # Reverse complement check (simplified - just look for palindromic sequences)
        rc_right = reverse_complement_simple(right)
        
        # Count matching positions
        matches = sum(1 for a, b in zip(left, rc_right) if a == b)
        inverted_score = matches / min(len(left), len(rc_right))
        
        # Combine scores (0-1 scale)
        combined_score = (gc_score * 0.4 + hpoly_score * 0.3 + inverted_score * 0.3)
        scores.append((i, combined_score))
    
    return scores


def reverse_complement_simple(seq: str) -> str:
    """Simple reverse complement for secondary structure analysis."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))


def identify_high_structure_regions(sequence: str, threshold: float = 0.7,
                                    window_size: int = 40) -> List[Tuple[int, int]]:
    """
    Identify regions likely to form stable secondary structures.
    
    Args:
        sequence: DNA sequence
        threshold: Score threshold for flagging regions (0-1)
        window_size: Window size for analysis
        
    Returns:
        List of (start, end) tuples for problematic regions
    """
    scores = simple_secondary_structure_score(sequence, window_size)
    
    regions = []
    in_region = False
    region_start = 0
    
    for pos, score in scores:
        if score >= threshold:
            if not in_region:
                region_start = pos
                in_region = True
        else:
            if in_region:
                regions.append((region_start, pos + window_size))
                in_region = False
    
    if in_region:
        regions.append((region_start, len(sequence)))
    
    return regions


def sequence_quality_score(sequence: str) -> dict:
    """
    Comprehensive quality scoring for a sequence.
    
    Returns dict with scores for various metrics (all 0-100 scale).
    """
    from .sequence_utils import gc_percent, longest_homopolymer
    
    sequence = sequence.upper()
    scores = {}
    
    # GC content (ideal 40-60%)
    gc = gc_percent(sequence)
    if 40 <= gc <= 60:
        scores["gc_content"] = 100
    elif 30 <= gc <= 70:
        scores["gc_content"] = 80
    else:
        scores["gc_content"] = max(0, 100 - abs(gc - 50))
    
    # Homopolymer runs (5bp is reasonable, >8bp is problematic)
    hpoly = longest_homopolymer(sequence)
    if hpoly <= 5:
        scores["homopolymer"] = 100
    elif hpoly <= 8:
        scores["homopolymer"] = 60
    else:
        scores["homopolymer"] = max(0, 100 - (hpoly - 8) * 10)
    
    # Repeat detection
    tandem = len(find_tandem_repeats(sequence))
    dispersed = len(find_dispersed_repeats(sequence))
    total_repeats = tandem + dispersed
    
    if total_repeats == 0:
        scores["repeats"] = 100
    elif total_repeats <= 3:
        scores["repeats"] = 80
    else:
        scores["repeats"] = max(0, 100 - total_repeats * 10)
    
    # Secondary structure risk
    high_struct = len(identify_high_structure_regions(sequence))
    if high_struct == 0:
        scores["secondary_structure"] = 100
    else:
        scores["secondary_structure"] = max(0, 100 - high_struct * 15)
    
    # Overall score (average of components)
    scores["overall"] = int(sum(scores.values()) / len(scores))
    
    return scores
