# Development log — decisions, discoveries, and bugs worth knowing about

This is a handoff document, not a changelog. It records the *reasoning* behind
non-obvious choices and the *root cause* of bugs that took real investigation
to track down — the kind of context that's expensive to rediscover but cheap
to read once. Organized by module/topic, not strictly chronologically.

## Phases 1-3 (mass engine, resolvability scoring, fragment ladder)

See `README.md`'s checklist for what's built. A few decisions worth flagging:

- **Mass calculation wraps `pyteomics` via `reticulate`** rather than
  reimplementing residue-mass tables — see `R/mass_calculation.R`. The Python
  venv is named `proteoformtracker-py` and is NOT project-relative (lives in
  reticulate's own managed location), so moving the project folder never
  breaks it.
- **PTM propagation onto relevant isoforms is NOT applied to the
  confounding-protein set** — PTM occupancy isn't a fixed proteome property
  the way sequence is; applying a target's PTM to unrelated confounders would
  overstate how distinguishable they are. See `ptracker_vibe.md`'s
  limitations section.
- **Confounding-protein search is a precomputed, offline mass index**
  (`R/reference_proteome_index.R`), never recomputed per query — cross-gene
  mass collisions are a proteome-scale lookup, not a per-target computation.
- **Middle-down confounders come from a protease-digested reference pool**,
  not the intact-protein index — a real collision risk for a 4 kDa target is
  another protein's 4 kDa *digest fragment*, not just another 4 kDa intact
  protein. First use of a given protease is a real ~20-40s proteome-wide
  digest; cached to disk after that (`get_digested_reference_pool()`).

## Option 3 (rMATS adapter) — the bulk of recent work

### The core problem this module solves

rMATS only reports a differential exon (or two, for MXE) and its immediate
flanking exons — never the rest of the transcript. There's no way to compute
a full protein from the event alone. The chosen strategy: look up which
*already-annotated* transcripts of the gene structurally match each "arm" of
the event (exon-inclusion vs. exon-skipping for SE; 1st-exon vs. 2nd-exon for
MXE) using the same precomputed exon index Option 1 uses, and hand the user
real, full-length transcripts rather than trying to invent one from a
handful of exons.

### rMATS' own coordinate conventions (confirmed empirically, not from docs)

- `*ES` columns are 0-based starts (need `+1`); `*EE` columns are already
  1-based inclusive.
- `upstream`/`downstream` column names are by **genomic coordinate** (lower/
  higher), not transcription direction. On a `-` strand gene this is the
  OPPOSITE of 5'/3' order. Confirmed directly against MYOM1 (a `-` strand
  gene): the "downstream" (higher-coordinate) flank has a LOWER exon_number
  (comes first, 5' side) than the target exon. See
  `genomic_flanks_in_transcript_order()`.

### Matching rule: why only ONE flank is required, not both

Original (wrong) design required both the arm's own exon's upstream AND
downstream neighbor to exactly equal rMATS' reported flanks. This produced
**zero matches** for real, useful cases — confirmed directly on IL32's MXE
event: the reported "1st exon" and its upstream flank had no exact match
anywhere in the gene's annotation, while the "2nd exon" and its downstream
flank matched 11 real transcripts (including the canonical one) exactly.
rMATS' own reported coordinates don't always land on an annotated boundary on
*both* sides — relaxing to "at least one flank must match AND be immediately
adjacent" recovers these real matches without weakening the evidence needed
(adjacency to even one flank is still strong, specific evidence that a
transcript represents that arm).

**This is evaluated independently per arm, per transcript, using that arm's
own exon** — not "does this transcript match somewhere." A transcript can
genuinely contain BOTH of MXE's "mutually exclusive" exons back-to-back
(confirmed directly: several real IL32 transcripts, e.g. `ENST00000878154`,
contain exon1 immediately followed by exon2). For such a transcript:
exon1's own neighbors are (unrelated exon, exon2) — neither is a flank, so it
fails the "1st exon form" test. exon2's own neighbors are (exon1, downstream
flank) — the downstream neighbor IS the flank, so it passes the "2nd exon
form" test. Both exons are really present in the same transcript; only one
arm's adjacency condition happens to be satisfied. This is *why* a highlight
box for exon1 can appear on a track that's only listed under the "2nd exon
form" arm — the visualization draws a box wherever a shown track's real,
full annotation contains that *exact* coordinate, independent of which arm
"claimed" that transcript as a candidate. Not a bug; two different, both
correct, definitions being shown side by side. (A UI note explaining this
now lives directly above the candidate list in `server.R`.)

### `reference_exon_index` is CDS-only — full exon structure needs a live fetch

The precomputed index is built from GTF `CDS` rows only (by design, for
Option 1's residue-numbering math). It has **no UTR exons at all** and can
clip a partially-coding exon's boundary to just its coding portion. Confirmed
directly: IL32's real, full first exon is `3065784-3065826`, but
`reference_exon_index`'s CDS-only row for the same locus is
`3065812-3065826` — missing the 5' UTR portion entirely, and this
*coincidentally* equals rMATS' own reported "1st exon" coordinate exactly
(rMATS reports against the full annotation, not the CDS-only subset).

Fix: the exon-alignment visualization fetches each real candidate's FULL
exon list via `fetch_transcript_exons()` (live Ensembl REST, disk-cached),
using `reference_exon_index` only as the CDS/coding overlay on top of that —
same pattern Option 2 already used for its own novel-vs-known comparison.

### Constructed ("synthetic") isoforms — building and translating a whole isoform

Every arm now gets a "Constructed" isoform track regardless of whether any
real transcript matched: splice a user-selectable **backbone** transcript's
own exons (everything OUTSIDE the local AS region, by genomic POSITION
thresholding — a backbone's own exon can legitimately extend further than
rMATS' reported flank boundary, e.g. real UTR fused onto a partially-coding
exon) together with rMATS' own reported local exons (flanks + cassette
exon(s)) taken verbatim. See `build_rmats_arm_isoform()`.

Default backbone preference order (`default_backbone_for_arm()`): (1) best
real match for this arm itself, (2) best-pairing-score real match from
ANOTHER arm of the same event (the common MXE case: one arm has real
matches, the other doesn't, but they're almost certainly variants of the
same underlying transcript), (3) gene's canonical transcript, (4) first
transcript of the gene. Always user-overridable via a per-arm dropdown.

**Translating the constructed isoform to a real protein** (added when the
user asked for constructed isoforms to be selectable for MS1/MS2, not just
shown in the visualization): extracts the genomic DNA for the isoform's own
CDS sub-ranges directly from the local genome FASTA (`samtools faidx
reference/genome/GRCh38.fa`, not a REST call — there's no Ensembl endpoint
for translating an arbitrary coordinate list), concatenates in ascending
genomic order, reverse-complements the WHOLE combined string if the
transcript is `-` strand (equivalent to, and simpler than, reverse-
complementing each fragment individually — `revcomp(A+B+C) ==
revcomp(C)+revcomp(B)+revcomp(A)`), then translates with a standard codon
table. **Validated by re-deriving a real transcript's own known Ensembl
protein sequence from its own CDS coordinates this way and diffing —
byte-for-byte identical on both a `+` and a `-` strand test transcript**
before wiring this into the UI.

One non-obvious translation detail: `reference_exon_index`'s CDS ranges
follow GTF/GENCODE convention and **exclude the terminal stop codon** — the
translator must NOT require finding an in-frame stop within the given
sequence (it won't be there for the common case); it translates the whole
thing and only truncates early if a stop genuinely appears mid-sequence.
Requiring a stop unconditionally was an early bug that silently discarded
every correctly-formed translation (see `translate_dna_standard()`'s doc
comment).

### Bugs found and fixed, with root causes

1. **Zero matches on real MXE data** — the both-flanks-required matching
   rule (see above). Fixed by relaxing to one flank + adjacency.
2. **"Missing exons" / single black box spanning all tracks** — initially
   misdiagnosed as a stale zoom/pan bug (real fix, kept, but not the actual
   cause). Real cause: CDS-only `reference_exon_index` used directly as the
   visualization's exon source, so UTR exons were simply never in the data.
   Fixed by switching real-candidate tracks to `fetch_transcript_exons()`
   (see above).
3. **SE mode's "Find matching transcripts" stayed permanently disabled**
   after switching from an MXE file to an SE file without clicking Reset.
   Root cause (confirmed via temporary server-side `message()` logging):
   Shiny's selectize widget's underlying VALUE literally does not change
   when the SAME string ("1") is set twice in a row across two unrelated
   event lists — the widget's *label* updates, but no `change` event fires,
   so the server never learns a new pick was made. Fixed by tagging the
   event dropdown's own option VALUES with a generation number (`"3_1"` not
   `"1"`) so re-picking "row 1" of a new list is always a genuinely
   different string. This let a separate, less robust "did the user pick
   something new" tracking flag be deleted entirely in favor of a direct
   string-prefix check against the current generation number.
4. **Picking a different backbone transcript appeared to do nothing** (kept
   silently reverting to the default) **and the constructed isoform's shown
   mass/AA count never updated** — these were the SAME bug, not two. Once
   `rmats_candidate_meta()`/`rmats_constructed()` (added to show AA/mass
   columns) started depending on the backbone dropdowns' own input values,
   `output$rmats_results_ui`'s `renderUI()` — which is what DRAWS those
   dropdowns — transitively depended on them too. So every backbone pick
   triggered a full re-render of the whole panel, which recreated the
   dropdown with `selected = default_backbone`, discarding the user's pick
   before the constructed isoform was ever rebuilt against it. Fixed by
   reading the dropdown's OWN prior value back via `isolate(input[[...]])`
   before falling back to the default — the same pattern the candidate
   checkboxes already used for their own `checked` state, just not yet
   applied to this widget. **General lesson: any `renderUI()` whose body
   reads a reactive that itself depends on an input the SAME `renderUI()`
   renders is a reset trap unless every such input's value is isolate()-read
   back as its own new default.**

## Known gaps / explicitly deferred

- Constructed isoforms feed MS1/MS2 fully now (translation + exon table +
  "Add to comparison" wiring all done), but **A3SS/A5SS/RI event types are
  not implemented** — only SE and MXE. The matching core
  (`match_rmats_arm_transcripts()` / `build_rmats_arms_result()`) was written
  generically enough that adding a new event type mainly means writing its
  own `parse_rmats_*()` + `match_rmats_*_transcripts()` pair, not touching
  the shared matching/visualization/translation code.
- A full reference-proteome-wide MS2 fragment collision index (as opposed to
  the current small explicit-set check) is scoped but not built.
