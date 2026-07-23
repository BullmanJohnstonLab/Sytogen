# SyToGen Decision Matrix — Field Reference

`decision_matrix.tsv` is the full audit trail of every edit SyToGen considered
while resolving restriction motifs. This document explains what each column
means, exactly how it's computed, and how to read the two different kinds of
rows it contains.

## Two kinds of rows

For each motif SyToGen tries to resolve, the decision matrix contains **one
of two things**:

1. **Candidate rows** — one row per valid edit SyToGen considered for that
   motif. If three synonymous codons could all silence a site, you get three
   rows, with `chosen = True` on exactly one of them (the winner) and
   `chosen = False` on the rest (documented alternatives, not failures).
2. **A single sentinel row** — if SyToGen could not construct *any* valid
   edit for a motif, you get exactly one row for it instead, with
   `skip_reason` populated and most of the per-candidate columns blank.

You can tell which kind of row you're looking at from `skip_reason`: blank
means a real candidate row, populated means a sentinel row.

## Column reference

| Column | Meaning | How it's derived |
|---|---|---|
| `motif` | The recognition sequence pattern (may include IUPAC ambiguity codes like `N`, `R`, `Y`) that this row is about. | Taken directly from the motif table you uploaded. |
| `motif_start` | 0-based start position of *this specific occurrence* of the motif in the sequence. | From motif discovery (`_parse_motifs`) — either the coordinates in your motif table, or found by searching the sequence for the pattern. |
| `motif_end` | 0-based end position of the occurrence, inclusive. | Same source as `motif_start`. Can be `>= sequence length` for a motif that spans the circular origin — that's intentional (see "Circular constructs" below), not an error. |
| `motif_strand` | `+` or `-` — which strand the occurrence was found on. | Recorded at discovery time based on which strand the pattern (or its reverse complement) matched. |
| `edit_position` | The 0-based genomic position of the single base this candidate changes. | The position of the one base that differs between the original codon/base and the replacement. Blank on sentinel rows (no edit was made). |
| `gene_id` | The gene this edit falls inside, if any. | Looked up via the position against annotated `CDS`/`ORF`/`Marker` features. Blank if the edit is in a non-coding region, or on sentinel rows. |
| `gene_strand` | `+`/`-` strand of that gene. | From the same feature lookup as `gene_id`. Blank under the same conditions. |
| `original_codon` | The original 3-base codon (coding regions) or single original base (non-coding regions) at the edit position. | Read directly from the sequence before editing. |
| `replacement_codon` | The proposed replacement codon or base. | The synonymous codon (coding) or alternate base (non-coding) this candidate proposes. |
| `AA_LetterCode` | One-letter amino acid code for the codon. | Standard genetic code translation of `original_codon`. Since every coding candidate is synonymous by construction, this is the same amino acid before and after. Blank for non-coding rows. |
| `synonymous` | `True` for every coding-region candidate row. | Coding candidates are only ever generated as synonymous single-base edits in the first place, so this is always `True` when it applies. Blank (not `False`) for non-coding rows, since the concept doesn't apply there — not because the edit failed to be synonymous. |
| `motifs_destroyed` | How many *specific, previously-existing* motif occurrences this exact edit destroys. | Checked per-occurrence: does the exact known span of each nearby registered motif still match its pattern (forward or reverse-complement) after the edit? Usually `1` (the motif this row is about) but can be higher if one edit happens to simultaneously destroy an overlapping/adjacent occurrence too. |
| `reasoning` | Plain-English explanation of this row. | Generated from whether the row was chosen, what changed, and (for coding rows) the codon-usage score that mattered. Deliberately doesn't restate the destroyed/created counts, since those already have their own columns — see them there instead of in the text. |
| `motifs_created` | How many *new* motif occurrences this edit introduces nearby. | **Always `0` for every row you'll actually see.** Any candidate that would create a new motif is marked invalid and discarded before it ever becomes a row — so a `motifs_created` value only shows up as the *reason* a candidate never appears, not as something you'll see printed as nonzero. The column exists for schema completeness/transparency about what was checked. If you want to know about newly introduced sites, see `new_motifs_check.json` (below) — a separate, whole-construct pass that catches anything the per-edit check might have missed. |
| `usage_score` | The codon-usage preference value for the replacement codon, from your codon usage table. | Looked up directly from the uploaded table (higher = more preferred). Always `0` for non-coding rows, since there's no codon-usage concept outside a gene. |
| `gc_preserving` | Whether this specific base substitution keeps the same GC-class as the original base (G↔C or A↔T, vs. crossing between them). | `True`/`False` computed directly from `original_codon`'s/`replacement_codon`'s differing base. Used only as a *tiebreaker* — see "How the winner is chosen" below. |
| `total_score` | The overall ranking score used to pick the winner among candidates for the same motif. | `(motifs_destroyed × 1000) + (usage_score × 100) − (1 × 10)`. The flat `−10` is a fixed per-edit cost (every candidate makes exactly one edit); it only matters when comparing across motifs, not between candidates for the same motif. |
| `chosen` | `True` on the one candidate actually applied to the sequence for this motif; `False` on every other valid alternative. | The candidate with the highest `total_score` wins; `gc_preserving` breaks an exact tie (see below). Exactly one `True` per motif that resolved successfully. |
| `skip_reason` | Populated only on sentinel rows — why no valid candidate existed at all. | One of: `blocked_by_protected_region`, `no_synonymous_codon`, `all_candidates_rejected`, `no_valid_edit`. See below for what each means. |
| `attempted_count` | On sentinel rows only: how many individual single-base edits were actually tried before concluding none worked. | Count of every substitution attempt made across all positions the motif spans. Blank on candidate rows. |
| `rejected_count` | On sentinel rows only: how many of those attempts were rejected. | Currently always equal to `attempted_count` — every attempt SyToGen counts here was, by definition of reaching a sentinel row, a rejected one. Kept as a separate column since it's conceptually distinct (and could diverge if the logic changes). Blank on candidate rows. |
| `top_rejection_reason` | On sentinel rows only: the single most common reason among all rejected attempts. | One of `Protected region`, `Not synonymous`, `Creates new motif`, or `does not destroy this motif` — see below. Blank on candidate rows. |
| `top_rejection_count` | On sentinel rows only: how many of the rejected attempts shared that top reason. | E.g. `12/12` worth of attempts all failing for the same reason is a strong signal, not a coincidence — see "Reading a sentinel row" below. Blank on candidate rows. |

## How the winner is chosen

For each motif with at least one valid candidate, SyToGen sorts candidates by:

1. **`total_score`**, highest first. This is dominated by `motifs_destroyed`
   (worth 1000 points each) — a candidate that destroys the motif always
   beats one that doesn't, and `usage_score` (worth up to ~100 points,
   depending on your codon table's scale) only matters as a second-order
   preference among candidates that are otherwise equally good at actually
   silencing the site.
2. **`gc_preserving`** breaks an exact tie in `total_score`. This never
   overrides a real difference in codon usage — it only decides between
   candidates that scored identically on everything else.

## Reading a sentinel (no-candidate) row

`skip_reason` tells you which of four situations applies:

- **`blocked_by_protected_region`** — every position the motif spans falls
  inside an annotated protected/regulatory feature. No edit was even
  attempted; there was nothing legal to try.
- **`no_synonymous_codon`** — the motif sits in a codon whose amino acid has
  no synonymous alternative (Met and Trp each have only one codon). Again,
  nothing was attempted.
- **`no_valid_edit`** — a catch-all for a position with no editable context
  at all (rare in practice).
- **`all_candidates_rejected`** — SyToGen *did* try edits (see
  `attempted_count`), and every single one was rejected. This is the one
  worth reading `top_rejection_reason`/`top_rejection_count` for:
  - If the top reason is **`does not destroy this motif`** and the count is
    100% of attempts (e.g. `12/12`), that's usually a sign of a *second
    copy of the same short motif* sitting nearby — the per-edit check looks
    in a window around the edit, and if a sibling occurrence of the same
    pattern is in that window, it can't tell "the target site is gone" from
    "a neighboring copy still matches." This isn't a dead end so much as a
    hint that this particular site may need a wider look.
  - **`Protected region`** — every attempted base fell inside a protected
    annotation.
  - **`Creates new motif`** — every attempted edit would introduce a new
    restriction site elsewhere nearby, so none were allowed through.
  - **`Not synonymous`** — the substitution would have changed the encoded
    amino acid.

## Circular constructs

For a motif that spans the circular origin (wraps from the end of the
sequence back to the start), `motif_end` is left un-clamped and can exceed
the sequence length — e.g. `motif_start: 2995, motif_end: 3004` on a
3000bp plasmid means the occurrence runs from position 2995 to the origin
and continues to position 4 (`3004 mod 3000`). This is the same convention
used internally for circular-aware position handling; it's not a bug, and
downstream tooling (assembly planning, primer design) already accounts for
it.

## Companion files

Two other output files reflect information adjacent to the decision matrix:

- **`summary.json`** — aggregate counts across the whole run (motifs found,
  resolved, unresolved, edits applied, and `new_motifs_introduced`).
- **`new_motifs_check.json`** — the result of the final, whole-construct
  safety check that runs once after all edits are applied: a full re-scan
  of the entire final sequence for every target pattern, independent of the
  per-edit checks that produced `motifs_created` above. If this reports
  anything nonzero, the SyToGen page will also show a warning banner after
  the run completes.
