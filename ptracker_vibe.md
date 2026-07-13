# ProteoformTracker — design spec

Planning tool for isoform/proteoform-level detectability in top-down proteomics (TDP). Companion project to **IsoPepTracker** (bottom-up). R Shiny.

**Reference implementation to fork/reuse:** https://github.com/HuangLabAtUAB/IsoPepTracker/

Reuse directly from IsoPepTracker rather than rebuilding:
- Input modules 1–3 (isoform selection informed by RNA-seq; long-read FASTA → ORF; rMATS event parsing into transcript/protein sequence)
- GTF/GFF parsing
- minimap2 alignment wrapper
- TransDecoder ORF-calling wrapper
- Exon-track visualization component (becomes the coordinate axis for the new fragment-ladder track, see Visualization outputs below)

## Scope

Answers a pairwise, comparative question — *given two specific proteoforms, can this instrument tell them apart* — not an absolute-size feasibility question. The tool does not compute or output a mass cutoff. Mass enters calculations only as it affects a specific pair's resolving power and charge-state envelope (see Orbitrap resolving power model).

Two things a naive Δmass-only check misses, which this tool must handle:
1. The Δmass needed for resolvability is not a fixed number — it scales with the pair's own mass and charge state via `M / R(m/z)`. Fixed-threshold rules ("Δmass ≥ 3 Da is enough") produce false negatives on large easily-separable pairs and false positives on small marginal pairs.
2. A pairwise-resolvable pair can still be confounded by an unrelated third protein (different gene) whose intact mass happens to fall within the resolving window of either member — the pair problem is never actually isolated from the rest of the proteome. This is the rationale for treating the MS2 fragment ladder as a first-class, independent axis of discrimination, not an optional extra: two proteins from different genes essentially never share extended sequence identity, so their fragment ions diverge quickly even when intact masses coincide.

## System overview

Three input modules → one shared `Proteoform` object → relevant-protein search + confounding-protein search → scoring engine → visualization outputs.

### The `Proteoform` object (shared internal schema)

All three input modules must resolve into this before any downstream code runs. Downstream modules (mass calculation, charge-envelope prediction, resolvability scoring, fragment-ladder generation) consume only this object, never raw module-specific input.

Minimum fields:
- Mature amino acid sequence (after N-terminal Met excision and signal peptide/propeptide removal)
- List of applied PTMs with Unimod-style mass deltas (empty by default)
- Provenance tag (which input module produced it)

### Relevant proteins vs. confounding proteins (implement as two separate searches)

- **Relevant proteins**: share exon structure/sequence similarity with the target. Same computation IsoPepTracker already does for exon-level comparison; can additionally use SQANTI3-style structural classification against a reference annotation for novel long-read/rMATS-derived isoforms.
- **Confounding proteins**: share similar intact mass regardless of sequence/exon relationship — candidates for cross-gene MS1/MS2 collisions. Implement as a fast, precomputed range lookup against a one-time offline calculation of theoretical intact masses for the full reference proteome. Do not recompute per user query.

## Technical modules

### 1. Proteoform mass calculation

Handle explicitly, beyond raw ORF translation:
- N-terminal Met excision (sequence-context dependent)
- Signal peptide/propeptide cleavage (mature form is what TDP measures)
- PTMs as mass deltas from a Unimod-style lookup table, populated only when user specifies known/suspected sites

Recommend wrapping an existing library (e.g. `pyteomics` via `reticulate`) rather than reimplementing mass tables.

### 2. Orbitrap resolving power model

Resolving power is quoted at one reference m/z (conventionally 200) and degrades as m/z increases:

```
R(m/z) = R_ref × √(m/z_ref / m/z)
```

Converts to a mass-domain peak width:

```
ΔM_FWHM (Da) ≈ M / R(m/z at the best available charge state)
```

**Worked example**: 16 kDa proteoform, z=20 → m/z ≈ (16,000 + 20×1.007)/20 ≈ 801. On R_ref=120,000 @ m/z_ref=200: R(801) = 120,000 × √(200/801) ≈ 60,000. ΔM_FWHM ≈ 16,000/60,000 ≈ 0.27 Da.

Apply a 1.5–2× safety margin over ΔM_FWHM before calling a pair confidently baseline-separable (not just theoretically resolved).

**Instrument-specific**: this scaling law is Orbitrap-only. FT-ICR follows a different relationship — do not reuse this formula if instrument-agnostic support is added later.

### 3. MS1 resolvability and envelope-crowding scoring

For a candidate pair (target vs. relevant isoform, or target vs. confounding protein):
1. Predict a plausible charge-state envelope for each proteoform (see Charge-state/S-N modeling).
2. Compute R(m/z) and ΔM_FWHM at every charge state in the envelope; take the charge state giving the smallest (best) ΔM_FWHM.
3. Compare best-case ΔM_FWHM × safety margin against the pair's actual Δmass → resolvable / marginal / not-resolvable verdict.
4. Separately check envelope crowding: plot each proteoform's predicted (m/z, charge) peak grid, flag regions where peaks from different proteoforms fall within one resolvable window of each other. Run against both the relevant-isoform set and the confounding-protein set.

Above ~25–30 kDa, isotope envelopes themselves widen enough that even a perfectly resolving instrument may not cleanly separate close proteoforms (envelopes interleave). This is an empirical threshold, not derived from the R(m/z) formula — validate against real example pairs (see Validation).

### 4. Charge-state envelope and signal-to-noise modeling

Condition both on a user-selectable acquisition mode (native vs. denatured):
- **Charge-state ceiling**: denaturing, ≤~40 kDa → tracks basic-residue count (sequence-derivable). Above that, or native-like generally → Rayleigh-limit relationship tied to surface area/mass.
- **S/N decay with mass**: total ion signal divides across more charge states and a wider isotope envelope as mass increases; per-peak S/N drops faster under denaturing than native conditions. Implement as a relative, normalized penalty against a user-supplied baseline — not an absolute intensity prediction (absolute signal depends on sample loading/instrument sensitivity the tool can't know).

### 5. Fragmentation propensity model

No trained intensity-prediction model (Prosit/MS2PIP-style) exists for intact-protein fragmentation — use a literature-grounded, rule-based propensity score instead. Four effects combine multiplicatively into a per-bond score:

- **Residue-pair effects**: cleavage N-terminal to proline and C-terminal to aspartate strongly enhanced under collision-based fragmentation; enhanced further under native vs. denatured top-down conditions.
- **Dissociation-method dependence**: residue-pair enhancements are method/mode-specific (majority of propensities differ significantly between native and denatured). ETD/ECD and UVPD produce more even, less residue-selective cleavage → flatter propensity table under those methods.
- **Positional (terminal-vs-internal) bias**: terminal fragment ions statistically favored to carry majority of fragment ion current for up to ~3 backbone cleavage events, regardless of protein size — weight cleavages near termini above deep-internal ones.
- **Structural accessibility (ETD/ECD only)**: efficiency depends on gas-phase charge density/higher-order structure; compact folded regions resist cleavage regardless of residue identity. Use a disorder/flexibility predictor (IUPred3) or AlphaFold DB per-residue confidence as a heuristic accessibility proxy — label clearly in UI as approximate, not a validated structural predictor.

**Key structural property**: because terminal fragment ladders are cumulative sums along the sequence, two isoforms' ladders are identical up to the sequence divergence point (typically a differential exon) and uniformly offset by Δmass beyond it — a cleaner pattern than bottom-up peptide uniqueness (which is scattered with no such structure).

### 6. Cross-protein MS2 fragment collision check

Two distinct reasons, flag separately in UI:
- **Database/scoring ambiguity**: a target fragment mass coinciding with a confounding protein's fragment is weaker confirmatory evidence, even without co-isolation.
- **Co-isolation/chimeric spectra**: if a confounding protein's intact mass is close enough to be co-isolated in the same precursor window, its fragments physically appear in the same MS2 scan (the scenario TopMPI-style tools disentangle post-acquisition). ProteoformTracker flags this risk pre-acquisition; it does not resolve chimeric spectra itself.

Because cross-gene proteins essentially never share extended sequence identity, most collisions should clear as resolvable by MS2 — this module's value is catching exceptions (short proteins, unusually homologous unrelated proteins, or discriminating fragments clustering in a low-propensity region).

**Implementation note**: reuses IsoPepTracker's peptide-uniqueness core logic (fragment/peptide mass unique against a defined background set), applied to terminal fragment ions instead of tryptic peptides, run against both the relevant-isoform set and the precomputed confounding-protein/reference-proteome fragment index.

## Visualization outputs

1. **Exon + fragment ladder track**: reuse IsoPepTracker's exon-track component as the shared coordinate axis; mirrored fragment-ladder rows (one per proteoform) replace the tryptic-peptide bars, tick opacity/color encodes propensity score, vertical marker at the sequence-divergence point.
2. **MS1 m/z envelope overlay**: interactive stick spectrum, predicted charge-state peaks for target/relevant isoforms/confounding proteins simultaneously, shaded band per peak = predicted FWHM width at chosen instrument resolution (overlap should be visually self-evident, not require reading numeric Δmass).
3. **Candidate feasibility scorecard**: one row per candidate pair — Δmass, MS1 resolvability, envelope-crowding flag, MS2 fragment-based discriminability, overall verdict.

## Third-party tools

Not repeated: GTF/GFF parsing, minimap2, TransDecoder, rMATS parsing (already in IsoPepTracker).

| Function | Tool(s) | Notes |
|---|---|---|
| Isoform structural classification ("relevant proteins") | SQANTI3; gffcompare/gffread | Classifies novel long-read/rMATS isoforms against reference annotation |
| Sequence similarity search | BLAST+ or DIAMOND | Fallback when exon-structure comparison unavailable |
| Mass calculation with PTMs | pyteomics (via reticulate); Unimod | N/C-terminal processing, mono/average mass, PTM deltas |
| Isotope envelope modeling | Averagine model; pyteomics.mass.isotopologues or IsoSpec | Averagine = field-standard approximation above a few kDa |
| Charge-state envelope simulation | UniDec (CLI/scriptable) | Forward-simulates charge-state distributions; upgrade path to real deconvolution validation |
| Fragment ladder generation | Custom (deterministic arithmetic), or pyteomics/pyOpenMS TheoreticalSpectrumGenerator | No external tool strictly required |
| Structural accessibility proxy (optional) | IUPred3; AlphaFold DB confidence | Heuristic only |
| Reference proteome/annotation | UniProt reference proteome; Ensembl/RefSeq GTF/GFF | Background for confounding-protein and fragment-collision searches |

## Implementation roadmap

**Phase 1 — core data model and mass engine**
- `Proteoform` schema
- Mass calculation with N/C-terminal processing + optional PTMs
- Offline precomputed reference-proteome mass index (confounding-protein lookup)

**Phase 2 — resolvability scoring and MS1 visualization**
- R(m/z) resolving-power model + ΔM_FWHM calculation
- Charge-state envelope prediction + S/N penalty model
- MS1 resolvability/envelope-crowding scorer + m/z envelope visualization

**Phase 3 — fragmentation and MS2 visualization**
- Fragment ladder generator + propensity scoring model
- Cross-protein fragment collision check
- Combined exon + fragment ladder track + candidate scorecard

**Phase 4 — validation**
- Benchmark set of real, published isoform pairs with known TDP outcomes → sanity-check scorer verdicts
- For borderline cases, generate synthetic spectra and run through UniDec to compare the closed-form envelope-crowding heuristic against real deconvolution behavior
- Document limitations (below) alongside validation results

## Known limitations (state plainly in tool docs, not implied solved)

- Envelope-crowding check (Section 3) is a closed-form heuristic, not a guarantee of real deconvolution-algorithm behavior — Phase 4 validates against UniDec.
- Fragmentation propensity model predicts cleavage likelihood, not observed ion intensity — no mature, generalizable intensity-prediction model exists for intact-protein fragmentation.
- S/N model is relative/normalized against a user-supplied baseline, not an absolute detectability prediction.
- Does not model gas-phase rescue strategies (proton transfer charge reduction, ion mobility) that can recover discrimination when standard MS1/MS2 falls short — a "marginal" verdict means "standard workflow likely insufficient," not "impossible under every workflow."
- Orbitrap R(m/z) scaling law is instrument-family specific — FT-ICR or other platforms need a different resolving-power relationship.
- PTMs are propagated from a target proteoform onto its relevant-isoform set (with site-position mapping via pairwise alignment, skipping isoforms where the site is absent or the mapped residue doesn't match — see `propagate_ptms_to_relevant_set()`), but never onto the confounding-protein set. PTM occupancy is sample/condition-specific biology, not a fixed proteome property, so there is no canonical "modified state" of an unrelated background protein to precompute — the confounding-protein search (both mass- and m/z-domain) always assumes unmodified reference sequences.
- Fragment ladder generation (Section 5) only models terminal ions (b/y, one backbone cleavage per fragment) — internal fragments (two simultaneous cleavages, anchored to neither terminus) are not modeled. This is a real phenomenon in top-down MS/MS, not fabricated, but is deliberately out of scope: the number of possible internal fragments is combinatorial in sequence length (`C(n-1,2)` vs. `2(n-1)` for the full terminal ladder — roughly 185x more for a 742-residue protein), and each individual internal fragment is inherently less probable than a terminal ion, since it requires two independent cleavage events rather than one — consistent with terminal ions dominating the fragment ion current (see Section 5's positional-bias effect).

## Future considerations (deferred, not yet built)

- **Confounder PTM caveat check**: a lightweight, opt-in sanity check that widens the confounding-protein search window by a user-specified common PTM mass delta (e.g. +80 Da for phospho) and reports whether that changes the candidate set — flags *whether* PTM-driven confounding could plausibly matter for a given target, without trying to enumerate which specific background proteins might carry the modification (which isn't knowable from sequence alone). Deferred in favor of documenting the limitation above; revisit if real validation cases show it's needed.

## Reference reading

- FLASHDeconv, TopFD, UniDec — open-source deconvolution engines
- TopPIC, ProSight, MSPathFinder — proteoform ID/database search engines
- TopMPI — identification of co-isolated, multiplexed proteoform spectra from two distinct proteins
- SQANTI3 — structural classification of long-read isoforms against reference annotation
- TransDecoder — ORF prediction from assembled/long-read transcripts
- Marty et al. — UniDec: Bayesian deconvolution of native/intact protein spectra
- Fornelli, Compton et al. — proton transfer charge reduction (PTCR) for overlapping precursor/product ion signals in top-down MS
- Studies on gas-phase fragmentation propensities under native/denatured top-down conditions (proline/aspartate cleavage enhancement)
- Studies on terminal vs. internal fragment ion production statistics in top-down MS
- Reviews of ion activation methods (HCD, ETD/ECD, UVPD) for top-down proteomics of large proteins

## Glossary — bottom-up (IsoPepTracker) → top-down (ProteoformTracker)

| Bottom-up concept | Top-down analogue |
|---|---|
| Tryptic peptide | Terminal fragment ion (b/y or c/z), generated in silico from the intact sequence, not proteolytic digestion |
| Peptide "isoform-unique" if no other isoform produces the same tryptic peptide | Fragment "isoform-specific" if it falls beyond the divergence point in the cumulative fragment ladder |
| Peptide mass/charge, typically 2+ to 3+ | Intact proteoform mass/charge, typically double- to triple-digit charge states |
| Peptide detectability is largely solved/static | Proteoform detectability is size- and instrument-dependent, must be actively modeled |
| Peptide-level database search | Proteoform-level search against intact + fragment masses (e.g. TopPIC, ProSight) |
| Digestion enzyme specificity (e.g. trypsin cleaves after K/R) | No enzymatic digestion; cleavage follows gas-phase fragmentation propensities |
| Peptide uniqueness checked against whole tryptic proteome | Fragment collision checked against confounding proteins + reference proteome |