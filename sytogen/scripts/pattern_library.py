"""
Pattern specification language and library for common biological sequences.

Provides named patterns for:
- Restriction sites
- Assembly standards (Golden Gate, BioBricks, etc.)
- Regulatory elements
- Synthetic biology markers
"""

import re
from typing import Dict, List, Optional, Set
from dataclasses import dataclass


@dataclass
class PatternSpec:
    """Specification for a named DNA pattern."""
    name: str
    pattern: str  # IUPAC regex or plain sequence
    description: str
    category: str  # "restriction", "assembly", "regulatory", "motif"
    is_regex: bool = False
    reverse_complement: bool = True  # Whether to also match RC
    
    def compile(self):
        """Compile pattern to regex."""
        if self.is_regex:
            return re.compile(self.pattern, re.IGNORECASE)
        else:
            # Convert IUPAC to regex
            iupac_regex = iupac_to_regex(self.pattern)
            return re.compile(iupac_regex, re.IGNORECASE)


class PatternLibrary:
    """Centralized library of named patterns."""
    
    def __init__(self):
        self.patterns: Dict[str, PatternSpec] = {}
        self._load_standard_patterns()
    
    def _load_standard_patterns(self):
        """Load built-in patterns."""
        
        # Common restriction enzyme sites
        restriction_sites = [
            ("EcoRI", "GAATTC", "EcoRI recognition site"),
            ("BamHI", "GGATCC", "BamHI recognition site"),
            ("HindIII", "AAGCTT", "HindIII recognition site"),
            ("XbaI", "TCTAGA", "XbaI recognition site"),
            ("SpeI", "ACTAGT", "SpeI recognition site"),
            ("PstI", "CTGCAG", "PstI recognition site"),
            ("SalI", "GTCGAC", "SalI recognition site"),
            ("KpnI", "GGTACC", "KpnI recognition site"),
        ]
        
        for name, seq, desc in restriction_sites:
            self.add_pattern(PatternSpec(
                name=name,
                pattern=seq,
                description=desc,
                category="restriction",
                reverse_complement=True
            ))
        
        # Golden Gate BsaI sites (for TypeIIS assembly)
        self.add_pattern(PatternSpec(
            name="BsaI_site",
            pattern="GGTCTC",
            description="BsaI recognition site (Golden Gate)",
            category="assembly",
            reverse_complement=True
        ))
        
        # BioBricks RFC10 standard parts
        self.add_pattern(PatternSpec(
            name="BioBricks_prefix",
            pattern="GAATTCGCGGCCGCTTCTAGAG",
            description="BioBricks RFC10 prefix (EcoRI-NotI-XbaI)",
            category="assembly"
        ))
        
        self.add_pattern(PatternSpec(
            name="BioBricks_suffix",
            pattern="TACTAGTAGCGGCCGCTGCAG",
            description="BioBricks RFC10 suffix (SpeI-NotI-PstI)",
            category="assembly"
        ))
        
        # Common regulatory elements
        self.add_pattern(PatternSpec(
            name="RBS_consensus",
            pattern="AGGAGG",
            description="Ribosome binding site (RBS) consensus",
            category="regulatory"
        ))
        
        self.add_pattern(PatternSpec(
            name="Stop_codon_TAA",
            pattern="TAA",
            description="Stop codon (TAA)",
            category="regulatory"
        ))
        
        self.add_pattern(PatternSpec(
            name="Stop_codon_TGA",
            pattern="TGA",
            description="Stop codon (TGA)",
            category="regulatory"
        ))
        
        self.add_pattern(PatternSpec(
            name="Stop_codon_TAG",
            pattern="TAG",
            description="Stop codon (TAG)",
            category="regulatory"
        ))
        
        # Problematic sequences to avoid
        self.add_pattern(PatternSpec(
            name="AAAAA",
            pattern="AAAAA",
            description="5+ adenine homopolymer (transcription termination)",
            category="motif"
        ))
        
        self.add_pattern(PatternSpec(
            name="TTTTT",
            pattern="TTTTT",
            description="5+ thymine homopolymer (transcription termination)",
            category="motif"
        ))
    
    def add_pattern(self, spec: PatternSpec):
        """Add a pattern to the library."""
        self.patterns[spec.name] = spec
    
    def get_pattern(self, name: str) -> Optional[PatternSpec]:
        """Get a pattern by name."""
        return self.patterns.get(name)
    
    def find_patterns(self, name: str, sequence: str) -> List[tuple]:
        """
        Find all occurrences of a named pattern in a sequence.
        
        Returns:
            List of (start, end, matched_sequence) tuples
        """
        spec = self.get_pattern(name)
        if not spec:
            raise ValueError(f"Unknown pattern: {name}")
        
        regex = spec.compile()
        matches = []
        
        for match in regex.finditer(sequence):
            matches.append((match.start(), match.end(), match.group()))
        
        # Also search reverse complement if enabled
        if spec.reverse_complement:
            rc_seq = reverse_complement_full(sequence)
            for match in regex.finditer(rc_seq):
                # Convert position back to forward strand
                fwd_start = len(sequence) - match.end()
                fwd_end = len(sequence) - match.start()
                matches.append((fwd_start, fwd_end, f"{match.group()}(RC)"))
        
        return sorted(matches)
    
    def pattern_names_by_category(self, category: str) -> List[str]:
        """Get all pattern names in a category."""
        return [name for name, spec in self.patterns.items() 
                if spec.category == category]
    
    def list_patterns(self) -> Dict[str, str]:
        """Return dict of pattern_name -> description."""
        return {name: spec.description for name, spec in self.patterns.items()}


def iupac_to_regex(iupac_pattern: str) -> str:
    """Convert IUPAC degenerate code to regex."""
    iupac_map = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
    }
    
    pattern = ""
    for base in iupac_pattern.upper():
        pattern += iupac_map.get(base, base)
    
    return pattern


def reverse_complement_full(sequence: str) -> str:
    """Full reverse complement."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement.get(base, base) for base in reversed(sequence.upper()))


# Global instance
_default_library = None


def get_default_library() -> PatternLibrary:
    """Get or create the default pattern library."""
    global _default_library
    if _default_library is None:
        _default_library = PatternLibrary()
    return _default_library
