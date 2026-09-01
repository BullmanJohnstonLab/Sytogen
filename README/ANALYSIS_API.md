# SyToGen Advanced Analysis API Documentation

## Overview

SyToGen now includes advanced sequence analysis features built on DNA Chisel principles. These endpoints provide tools for detecting problematic sequences, analyzing quality metrics, and finding known patterns in DNA sequences.

## New Endpoints

### Sequence Analysis Endpoints

#### 1. `/api/analyze/repeats` - Detect Tandem and Dispersed Repeats

**Method:** POST

**Description:** Analyzes a sequence for both tandem repeats (consecutive identical motifs) and dispersed repeats (same sequence at multiple non-adjacent locations).

**Request:**
```json
{
  "sequence": "ATATATGCGCGCATATATAT",
  "min_repeat_length": 4  // optional, default 4
}
```

**Response:**
```json
{
  "tandem_repeats": [
    {
      "sequence": "AT",
      "positions": [0, 2, 4, 6, 14, 16, 18],
      "copy_count": 4,
      "length": 2
    }
  ],
  "dispersed_repeats": [
    {
      "sequence": "ATAT",
      "positions": [0, 14],
      "copy_count": 2
    }
  ]
}
```

**Use Cases:**
- Identify recombination hotspots (dispersed repeats)
- Detect homopolymer runs that might cause slippage
- Find sequences prone to tandem duplication events

**Limits:** Max 50,000 bp

---

#### 2. `/api/analyze/secondary-structure` - Identify High-Risk Secondary Structure Regions

**Method:** POST

**Description:** Analyzes sequence for regions prone to stable secondary structures (stem-loops, hairpins). Uses heuristic scoring based on GC content, homopolymer presence, and palindromic sequences.

**Request:**
```json
{
  "sequence": "GCGCGCGCGCGCATATATAT",
  "threshold": 0.7  // optional, 0-1 scale, default 0.7
}
```

**Response:**
```json
{
  "high_risk_regions": [
    {
      "start": 0,
      "end": 12,
      "length": 12
    },
    {
      "start": 18,
      "end": 28,
      "length": 10
    }
  ],
  "threshold": 0.7,
  "total_risky_bp": 22
}
```

**Scoring Heuristics:**
- GC content: Peak stability at 50% (0% and 100% GC both score low)
- Homopolymer runs: Long runs (10+ bp same base) increase risk
- Inverted repeats: Palindromic sequences prone to hairpin formation
- Returns stability score 0-1 per position

**Use Cases:**
- Avoid regions likely to misfold during synthesis
- Identify hairpin formation sites
- Optimize sequences for expression without secondary structure interference

**Limits:** Max 50,000 bp

---

#### 3. `/api/analyze/quality` - Comprehensive Sequence Quality Score

**Method:** POST

**Description:** Generates a complete quality assessment of a sequence based on 5 metrics.

**Request:**
```json
{
  "sequence": "ATGCATGCATGC"
}
```

**Response:**
```json
{
  "gc_content": 100,
  "homopolymer": 100,
  "repeats": 100,
  "secondary_structure": 85,
  "longest_homopolymer_length": 1,
  "overall": 93
}
```

**Metrics (0-100 scale):**
- **gc_content**: Scoring function peaking at 40-60% GC. 100% if in range, otherwise scaled.
- **homopolymer**: 100 if max run ≤5bp, 60 if ≤8bp, scales below that
- **repeats**: 100 if no repeats, 80 if 1-3 repeats, scales with more repeats
- **secondary_structure**: 100 if no high-risk regions, scales with detected structure risk
- **overall**: Simple average of all 4 component scores

**Use Cases:**
- Pre-screening sequences before synthesis
- Comparing quality of design alternatives
- Setting acceptance criteria for generated sequences

**Limits:** Max 50,000 bp

---

### Pattern Library Endpoints

#### 4. `/api/patterns/list` - List All Available Patterns

**Method:** GET

**Description:** Returns all patterns in the library organized by category.

**Response:**
```json
{
  "restriction": {
    "EcoRI": "EcoRI recognition site",
    "BamHI": "BamHI recognition site",
    "HindIII": "HindIII recognition site",
    "XbaI": "XbaI recognition site",
    "SpeI": "SpeI recognition site",
    "PstI": "PstI recognition site",
    "SalI": "SalI recognition site",
    "KpnI": "KpnI recognition site"
  },
  "assembly": {
    "BsaI_site": "BsaI recognition site (Golden Gate)",
    "BioBricks_prefix": "BioBricks RFC10 prefix (EcoRI-NotI-XbaI)",
    "BioBricks_suffix": "BioBricks RFC10 suffix (SpeI-NotI-PstI)"
  },
  "regulatory": {
    "RBS_consensus": "Ribosome binding site (RBS) consensus",
    "Stop_codon_TAA": "Stop codon (TAA)",
    "Stop_codon_TGA": "Stop codon (TGA)",
    "Stop_codon_TAG": "Stop codon (TAG)"
  },
  "motif": {
    "AAAAA": "5+ adenine homopolymer (transcription termination)",
    "TTTTT": "5+ thymine homopolymer (transcription termination)"
  }
}
```

**Categories:**
- **restriction**: Restriction enzyme recognition sites (8 sites)
- **assembly**: Assembly standard sites and prefixes (3 patterns)
- **regulatory**: Regulatory elements and stop codons (4 patterns)
- **motif**: Problematic sequences to avoid (2 patterns)

**Total:** 17 built-in patterns

---

#### 5. `/api/patterns/search` - Search Named Patterns in Sequence

**Method:** POST

**Description:** Finds all occurrences of named patterns in a sequence. Automatically searches both forward and reverse complement strands (when applicable).

**Request:**
```json
{
  "sequence": "ATGAATTCGGATCCGC",
  "patterns": ["EcoRI", "BamHI", "RBS_consensus"]
}
```

**Response:**
```json
{
  "EcoRI": [
    {
      "start": 3,
      "end": 9,
      "matched": "GAATTC"
    }
  ],
  "BamHI": [
    {
      "start": 9,
      "end": 15,
      "matched": "GGATCC"
    }
  ],
  "RBS_consensus": []
}
```

**Features:**
- Reverse complement matching for palindromic sites
- IUPAC degenerate code support (N, R, Y, S, W, K, M, B, D, H, V)
- Multiple patterns in single request
- Position-based matching with 0-based indexing

**Use Cases:**
- Finding restriction sites in constructs
- Validating assembly standards in design
- Detecting regulatory elements
- Scanning for problematic sequences

**Errors:**
- Returns `{"pattern_name": {"error": "Unknown pattern: ..."}}` for unknown patterns

---

## Integration with SyToGen Pipeline

### Workflow

1. **Before Synthesis Design**: Use `/analyze/quality` to pre-screen constructs
2. **During Design Iteration**: Use `/patterns/search` to detect unwanted sites
3. **Optimization Step**: Use `/analyze/repeats` and `/analyze/secondary-structure` to identify problematic regions
4. **Report Generation**: Use endpoints with internal report generation for documentation

### Example Workflow

```python
# 1. Upload construct
construct = "ATGCNNNATGC..."

# 2. Check overall quality
quality = POST /analyze/quality with construct
if quality["overall"] < 75:
    # Flag for manual review

# 3. Scan for restriction sites
sites = POST /patterns/search with ["EcoRI", "BamHI", ...]
if unexpected_sites:
    # Redesign to avoid

# 4. Detect problematic repeats
repeats = POST /analyze/repeats
if repeats["dispersed_repeats"]:
    # Flag recombination hotspots

# 5. Identify secondary structures
structures = POST /analyze/secondary-structure
if structures["total_risky_bp"] > 500:
    # Redesign high-risk regions
```

---

## Error Handling

All endpoints return standard error responses:

**400 Bad Request:**
```json
{
  "error": "Sequence required"
}
```

**413 Payload Too Large:**
```json
{
  "error": "Sequence too long (max 50,000 bp)"
}
```

**500 Server Error:**
```json
{
  "error": "Error message"
}
```

---

## Rate Limiting

These endpoints use **LIGHT_RATE_LIMIT** (1000 requests/hour per IP).

---

## Future Enhancements

### Planned Features

1. **Codon Synthesis Problems** (/api/optimize/codon-synthesis)
   - Template-based sequence optimization
   - Constraint satisfaction solving
   - Multi-objective optimization (GC%, codon usage, motif avoidance)

2. **Advanced Report Generation** (/api/report/<job_id>)
   - HTML reports with before/after comparisons
   - Interactive visualizations
   - PDF export capability

3. **Pattern Composition Language**
   - Combine multiple patterns with logical operators
   - Save custom patterns
   - Named pattern templates

### Extensibility

- Add custom patterns via library registration
- Extend scoring heuristics
- Integrate with DNA Chisel optimization engine

---

## Implementation Details

### Modules

- **dna_analysis.py**: Core repeat detection, secondary structure scoring
- **pattern_library.py**: Pattern management and matching
- **report_generator.py**: Report formatting and export

### Performance

- Repeat detection: O(n²) for tandem, O(n*m) for dispersed (n=sequence length, m=motif length)
- Pattern matching: O(n) per pattern via regex
- Quality scoring: O(n) with 4-bp window scanning
- Max input: 50,000 bp

### Dependencies

- BioPython 1.88 (sequence manipulation)
- Regex module (pattern matching)
- Plotly 6.9.0 (report visualization, optional)

---

## Testing

Complete test suite in `tests/test_analysis_modules.py`:
- 31 comprehensive tests
- Coverage: repeat detection, secondary structure, quality scoring, pattern library, report generation
- All tests passing ✅

---

## Example Use Cases

### Case 1: Quality Control Before Synthesis

```json
POST /analyze/quality
{
  "sequence": "ATGATGATGATGATGATG..."
}
```

Use to verify construct meets minimum quality threshold before ordering synthesis.

### Case 2: Identifying Restriction Sites in Design

```json
POST /patterns/search
{
  "sequence": "ATGCAATTCGGATCCATGC...",
  "patterns": ["EcoRI", "BamHI", "HindIII", "XbaI"]
}
```

Verify that planned restriction sites are present, unwanted sites are absent.

### Case 3: Detecting Problematic Repeats

```json
POST /analyze/repeats
{
  "sequence": "...construct sequence..."
}
```

Identify recombination hotspots (dispersed repeats) that could cause instability.

### Case 4: Secondary Structure Analysis

```json
POST /analyze/secondary-structure
{
  "sequence": "...construct sequence...",
  "threshold": 0.6
}
```

Find regions prone to hairpin formation that might interfere with expression.
