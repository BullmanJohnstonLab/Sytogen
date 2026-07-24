#!/usr/bin/env python3

# Import statements are grouped by standard library, third-party, and project-specific imports.
# Standard library
import copy
import math
import re

# Third-party
import Bio
from Bio.Seq import Seq

# Project constants
from ..constants import (
    DEFAULT_MONOVALENT_CONC_M,
    DEFAULT_OLIGO_CONC_M,
    GAS_CONSTANT_CAL,
    NN_INIT_AT,
    NN_INIT_GC,
    NN_THERMO_PARAMS,
)

# =============================================================================
# Basic sequence operations
# =============================================================================

# Basic sequence manipulation functions, including reverse complement, translation, and GC content calculation.

_RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    if not isinstance(seq, str):
        raise TypeError("Sequence must be provided as a string.")
    if any(base not in "ACGTacgt" for base in seq):
        raise ValueError("Sequence contains non-canonical bases.")
    return seq.translate(_RC_TABLE)[::-1]

def translate_sequence(seq: str) -> str:
    return str(Seq(seq).translate())

def translate_codons(input_sequence, codon_usage_table):
    translated_sequence = ''.join([str(codon_usage_table[x]['Translation']) for x in re.split(
        r'(\w{3})', input_sequence) if x != ''])
    return translated_sequence

def extract_subsequence(sequence, irange_):
    '''
    irange must be output form irange function
    '''
    # irange_ is expected to be an iterable of integer indices
    return "".join(sequence[i] for i in irange_)

def extract_window(sequence: str, start: int, length: int, circular: bool) -> str:
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

def extract_window_reverse(sequence: str, end_exclusive: int, length: int, circular: bool) -> str:
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

# =============================================================================
# Composition
# =============================================================================

def gc_percent(seq: str) -> float:
    seq = seq.upper()

    if not seq:
        return 0.0

    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)

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
# Thermodynamics
# =============================================================================

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

# =============================================================================
# Pattern searching
# =============================================================================

def find_patterns(
    sequence: str,
    patterns: list[str],
    circular: bool = True,
) -> dict[str, list[tuple]]:
    '''
    This function takes as input the sequence of the plasmid, and the list of target motifs (a list of strings)
    and returns a dictionary with keys each target motifs.
    For each key in the output a list of the following elements is provided:

    -the sequence on the plasmid corresponding to the target motif
    -the starting position of the target (included)
    -the ending position of the target (excluded)
    -the strand on which it is found (1 indicates the forward, -1 indicates the reverse complement)
    -junction/inner if the sequence is found on the 0 position or inside the actual sequence of the plasmid
    '''

    def plan_ambiguity(pattern: str) -> str:
        val = Bio.Data.IUPACData.ambiguous_dna_values
        re_pattern = ''
        for el in pattern:
            re_pattern = re_pattern + '[' + val[el] + ']'
        return re_pattern

    patts: dict[str, list[tuple]] = {}

    for pattern in patterns:
        tmp_fw = plan_ambiguity(pattern)
        data = sequence + sequence[:len(pattern)] if circular else sequence

        patts[pattern] = [
            (
                data[j.start() % len(sequence):j.end() % len(sequence)],
                j.start() % len(sequence),
                j.end() % len(sequence),
                1,
            )
            for j in re.finditer(tmp_fw, data)
        ]

        rc_sequence = str(Seq(data).reverse_complement())
        patts[pattern] += [
            (
                data[(len(data) - j.end()) % len(sequence):(len(data) - j.start()) % len(sequence)],
                (len(data) - j.end()) % len(sequence),
                (len(data) - j.start()) % len(sequence),
                -1,
            )
            for j in re.finditer(tmp_fw, rc_sequence)
        ]

    return patts

# =============================================================================
# Circular range handling
# =============================================================================

def circular_range(i, j, by=1, l=0):
    '''

    This function computes the range for a circular sequence.
    This function takes as input three integers:

    i: starting position of the range (included)
    j: ending position of the range (excluded)
    l: length of the sequence. The starting position of the sequence is supposed to be 0.

    The assumption is that the ranges are indicated according to c/python language convention:
    - numbering starts with 0
    - last position of a sequence is the actual length of the sequence reduced by one
    - ranges are made with first position included, second excluded

    The output is a list of lists with:
    - a single element (list) that corresponds to the actual range if the first position is lower than the second
    - two elements (lists) if the first position is higher than the second, indicating that we are considering
    the section of the circular sequence that goes from before the end of the sequence, across the junction
    to the starting portion of the sequence

    Example:

    >>> ss = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    >>> intervals = [[10, 1], [0, 2], [8, 9], [10, 11]]

    >>> intervals = [irange(10, 2, len(ss)), irange(-2, 3, len(ss)), irange(5, 9, len(ss)), irange(8, 11, len(ss))]
    >>> info(intervals)
    [[[10, 11, 12, 13, 14, 15], [0, 1]], [[14, 15], [0, 1, 2]], [[5, 6, 7, 8]], [[8, 9, 10]]]

    '''

    i_ = i % l
    j_ = j % l

    if i_ <= j_:
        return list(range(i_, j_, by))
    else:
        return list(range(i_, l, by)) + list(range(0, j_, by))

def merge_overlapping_ranges(list_ranges_):
    '''
    This function defines overlapping ranges and return a list of non overlapping ranges, in which the union of those overlapping is considered.
    In SyToGen, each irange must be elongate on the left and on the right by the length of the longest target motif. This must be done when defining
    the inputs of the irange function irange(i-m, j+m, l). This permits to retrieve from this function slices' positions that are at most contiguous
    and non-overlapping. Contigous slices preserve the correctness of the slice evaluation because of the two tails introduced by the irange.

    In this respect this function permits to define the positions of each bases of the slice on the plasmid.
    To get the irange take the first position and the last + 1 for each list in the list.

    >>> for el in list:
    >>>     info( irange(el[0][0], el[-1][-1]+1)) <-- This is the correct irange

    Each slice is made by:

    ---n [[], [], []] n--- :
    -two external n tails from the original sequence, where n is the length of the longest target motif
    -a list of lists of inner positions (from handle_ranges)

    The input is a list of iranges.

    Example:

    >>> ss = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    >>> intervals = [[10, 1], [0, 2], [8, 9], [10, 11]]

    >>> intervals = [irange(10, 2, len(ss)), irange(-2, 3, len(ss)), irange(5, 9, len(ss)), irange(8, 11, len(ss))]
    >>> info(intervals)
    [[[10, 11, 12, 13, 14, 15], [0, 1]], [[14, 15], [0, 1, 2]], [[5, 6, 7, 8]], [[8, 9, 10]]]

    >>> out_intervals = handle_ranges(intervals)

    >>> info(out_intervals)
    [[[8, 9, 10, 11, 12, 13, 14, 15], [0, 1, 2]], [[5, 6, 7, 8]]]

    '''

    list_ranges = copy.deepcopy(list_ranges_)
    set_tot_overlap = [list_ranges[0]]

    for i in range(1, len(list_ranges)):
        overlap = False
        for j in range(len(set_tot_overlap)):

            for k in range(len(list_ranges[i])):
                for l in range(len(list_ranges[j])):
                    if set.intersection(set(list_ranges[i][k]), set(list_ranges[j][l])):
                        set_tot_overlap[j][l] = list(
                            set.union(set(list_ranges[i][k]), set(list_ranges[j][l])))
                        overlap = True

        if not overlap:
            set_tot_overlap.append(list_ranges[i])

    return set_tot_overlap

# =============================================================================
# IUPAC pattern compilation
# =============================================================================

IUPAC_MAP = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]"}

def compile_iupac(motif):
    pattern: str = "".join(IUPAC_MAP[b] for b in motif.upper())
    return re.compile(pattern)


GC_BASES = frozenset("GC")
AT_BASES = frozenset("AT")


def is_gc_preserving_swap(old_base, new_base):
    """
    True if old_base -> new_base stays within the same GC-class (a G<->C
    or A<->T swap) rather than crossing between them. A same-class swap
    never changes local GC content at all.

    This ranks BELOW codon-usage preference (see run_sytogen_pipeline's
    candidate sort in sytogen_runner.py) — it's a tiebreaker for candidates
    that are otherwise equally good, not a reason to pick a worse-usage
    codon over a better one.
    """
    if not isinstance(old_base, str) or not isinstance(new_base, str):
        return False

    old_base = old_base.upper()
    new_base = new_base.upper()

    if len(old_base) != 1 or len(new_base) != 1:
        return False

    if old_base not in "ACGT" or new_base not in "ACGT":
        return False

    return (
        (old_base in GC_BASES and new_base in GC_BASES)
        or (old_base in AT_BASES and new_base in AT_BASES)
    )

