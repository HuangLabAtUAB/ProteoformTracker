# ProteoformTracker

An interactive online tool for planning proteoform detectability in top-down and middle-down
proteomics (TDP/MDP). Given a set of proteoforms — a specific combination of protein isoform and
post-translational modifications (PTMs) — ProteoformTracker prospectively models whether a mass
spectrometer can actually tell them apart, before any instrument time is committed.

ProteoformTracker is freely available at [https://www.proteoformtracker.org/](https://www.proteoformtracker.org/)

Tutorial & Documents: [https://huanglabatuab.github.io/proteoformtracker-docs/](https://huanglabatuab.github.io/proteoformtracker-docs/)

Companion to [IsoPepTracker](https://github.com/HuangLabAtUAB/IsoPepTracker/) (bottom-up).

## What it does

Two related questions, both answered before acquisition:

1. **Relevant-proteoform comparison.** Given a set of proteoforms you already care about (e.g.
   isoforms of the same gene, with or without PTMs), how separable are they from each other in
   MS1 (charge-state envelope, isotope pattern) and MS2 (b/y fragment ladder)?
2. **Confounding-protein search.** For a chosen target proteoform, does any *unrelated* protein in
   the full human proteome share its intact mass or a charge-state m/z peak closely enough to be a
   real confounder — a risk invisible if you only ever compare your proteoforms against each
   other?

### Input paths

Any of three paths resolves into a shared internal `Proteoform` representation (sequence + PTMs)
before any scoring runs:

- **Gene → isoform** — a gene's annotated Ensembl isoforms, fetched live via Ensembl REST
  (disk-cached).
- **FASTA → ORF** — a pasted/uploaded transcript sequence, spliced-aligned against GRCh38 with
  minimap2 to identify the matching known isoform, independently translated via TransDecoder ORF
  calling.
- **rMATS event** — an uploaded rMATS SE/MXE/RI/A5SS/A3SS differential-splicing results file; each arm of a
  selected event is matched against a precomputed genome-wide exon-structure index to recover
  full-length transcripts (including constructed synthetic isoforms for the alternative arm).

PTMs are specified per proteoform through a validated text syntax
(`<residue position>_<amino acid>_<Unimod name>`) drawn from a curated Unimod lookup table.

### MS1: charge-state envelope and isotope pattern

Predicts a plausible charge-state envelope (basic-residue-count ceiling for small denatured
proteoforms; Rayleigh-limit scaling above ~40 kDa or under native conditions), then computes the
**exact** isotope pattern from the proteoform's real elemental composition (via `pyteomics` +
IsoSpecPy), not a mass-only Averagine approximation. Predicted isotope spacing and pairwise
peak collisions are evaluated against an Orbitrap resolving-power model,
R(m/z) = R<sub>ref</sub> × √(m/z<sub>ref</sub>/m/z), with a user-adjustable safety margin.

### MS2: fragment ladder and fragmentation propensity

Generates the full deterministic b/y fragment-ion ladder and scores each backbone bond's
fragmentation propensity two ways:

- **Calibrated (default)** — a multiplicative formula (residue-pair, terminal-position,
  local-charge-density, proteoform-length, and phospho-proximity terms) calibrated by pooled
  logistic regression against ~2.7 million matched b/y-ion observations from two independent
  public top-down datasets (MassIVE MSV000094311, MSV000098558).
- **Length-free RF ranking** — a random-forest model trained on the same calibration data,
  excluding the length term, for within-proteoform relative ranking on very long proteoforms where
  the calibrated formula's length penalty can suppress every bond below fixed tier thresholds.

### Confounding-protein search

Searches a precomputed offline index of the full reference human proteome (UniProt reviewed
canonical, *n* = 20,391) on two independent axes: a mass-domain window search and an m/z-domain
charge-state collision search (catches proteins whose intact mass is far from the target but whose
charge-state peaks still collide — e.g. two masses in a 4:3 ratio colliding at paired charge
states). Because the m/z-domain search alone can return hundreds of real hits, the search step
stays cheap (no isotope computation) and returns a capped candidate list with a lightweight
pairwise MS2 shared-ion count; only after the user curates which candidates are worth
including — using outside evidence such as RNA-seq expression — does the full isotope-level
comparison run.

### Middle-down support

Simulates limited, partial in-silico digestion (OmpT, Lys-C, Lys-N, Glu-C, Asp-N) of every checked
proteoform within a user-set mass window, ranks candidate peptides by feasibility (missed-cleavage
count, predicted MS1 peak width) and PTM-site coverage, and applies the identical MS1/MS2 modeling
and confounder search to whichever candidates are selected — searched against other proteins' own
digest peptides, not intact proteins.

### Visualization

Results render as interactive, client-side SVG charts (`www/ptracker_viz.js`, vanilla JavaScript,
no charting library) with pan/zoom/hover/click, live resolvability statistics, and exon-coordinate
alignment across compared proteoforms.

##

Online portal is recommended: [https://www.proteoformtracker.org/](https://www.proteoformtracker.org/)

## Setup (not necessary if online portal is available)

Requires R ≥ 4.0 and Python (a project-local virtualenv is created automatically).

```r
# 1. Restore R package versions
renv::restore()

# 2. rtracklayer (Bioconductor) isn't in renv.lock -- only needed to build the exon index
BiocManager::install("rtracklayer")
```

```sh
# 3. Create the Python venv + install pinned pyteomics (<5.0, see python/requirements.txt for why)
Rscript scripts/setup_python_env.R
```

External CLI tools (via Homebrew), needed for the FASTA and rMATS input paths:
`samtools`, `minimap2`, `TransDecoder.LongOrfs` / `TransDecoder.Predict`.

Reference data (gitignored -- regenerate or copy separately; see `CLAUDE.md` for the full table
and download commands):

```sh
# Reference human proteome + mass index
Rscript scripts/build_reference_proteome.R

# m/z-domain charge-state collision index
Rscript scripts/build_mz_collision_index.R

# Genome-wide exon structure index (requires an Ensembl GTF + rtracklayer)
Rscript scripts/build_exon_index.R
```

Run the test suite:

```sh
Rscript tests/testthat.R
```

Run the app:

```r
shiny::runApp(port = 7412, host = "127.0.0.1")
```

## Project layout

- `ui.R` / `server.R` / `global.R` — the Shiny app; `global.R` sources every `R/*.R` file and loads
  precomputed indexes at startup.
- `R/proteoform_schema.R` — the shared `Proteoform` object every input path resolves into.
- `R/gene_isoform_data.R`, `R/ensembl_protein_fetch.R` — gene/isoform input path.
- `R/fasta_pipeline.R` — FASTA input path (minimap2 + TransDecoder).
- `R/rmats_adapter.R` — rMATS event input path.
- `R/mass_calculation.R`, `R/isotope_envelope.R`, `R/charge_envelope.R`, `R/ms1_scoring.R`,
  `R/resolving_power.R` — MS1 modeling.
- `R/fragment_ladder.R`, `R/fragmentation_propensity.R`, `R/fragmentation_propensity_rf.R`,
  `R/fragment_collision.R` — MS2 modeling.
- `R/reference_proteome_index.R`, `R/mz_collision_index.R` — confounder search indices.
- `R/digestion.R` — middle-down protease digestion.
- `R/viz_json.R`, `R/svg_render.R`, `www/ptracker_viz.js` — visualization payload assembly and
  client-side rendering.
- `tests/testthat/` — run via `Rscript tests/testthat.R`.

## License

MIT — see [LICENSE](LICENSE).
