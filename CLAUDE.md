# ProteoformTracker — working notes for Claude Code

R Shiny app that answers a pairwise question for top-down proteomics: *given two
specific proteoforms, can this instrument tell them apart?* Companion to
[IsoPepTracker](https://github.com/HuangLabAtUAB/IsoPepTracker/) (bottom-up).

Read these first, in this order:
1. `ptracker_vibe.md` — the original design spec (scope, scoring model, why the
   tool is pairwise/comparative rather than an absolute mass-cutoff check).
2. `README.md` — phase-by-phase feature checklist (what's done vs. planned).
3. `docs/DEVLOG.md` — narrative log of non-obvious decisions and bugs fixed
   during development, with the reasoning behind them. Read this before
   touching `R/rmats_adapter.R`, the exon-alignment renderer, or anything
   involving Shiny `renderUI()` re-render loops — several real bugs there
   only make sense with that context.

## Running it

```r
shiny::runApp(port = 7412, host = "127.0.0.1")
```
(`.claude/launch.json` already wires this up for Claude Code's browser-preview
tooling under the name `shiny-app`.)

## Dependencies (one-time setup)

- **R packages**: `renv::restore()` — restores everything in `renv.lock`.
- **`rtracklayer`** (Bioconductor): NOT in `renv.lock` (installed via
  `BiocManager::install("rtracklayer")` outside the renv snapshot — a real gap,
  not intentional). Only needed to build the exon index
  (`scripts/build_exon_index.R`); the running app doesn't import it directly.
  Install manually: `BiocManager::install("rtracklayer")`.
- **Python / reticulate**: `Rscript scripts/setup_python_env.R` creates the
  `proteoformtracker-py` virtualenv and installs `python/requirements.txt`
  (pyteomics, pinned `<5.0` — see that file for why). Mass calculation
  (`R/mass_calculation.R`) wraps `pyteomics` via `reticulate` rather than
  reimplementing residue-mass tables.
- **External CLI tools** (via Homebrew): `samtools`, `minimap2`,
  `TransDecoder.LongOrfs` / `TransDecoder.Predict`. Required for Option 2
  (FASTA→ORF alignment, `R/fasta_pipeline.R`) and Option 3's constructed-isoform
  translation (`R/rmats_adapter.R`, uses `samtools faidx` directly on the
  genome FASTA — see DEVLOG for why this exists instead of a REST call).

## Reference data (gitignored — `/data/*` and `/reference/` — must be
regenerated or copied separately; they don't travel with `git clone`)

| Path | How to get it |
|---|---|
| `data/annotation/Homo_sapiens.GRCh38.*.gtf.gz` | `curl -o data/annotation/Homo_sapiens.GRCh38.116.gtf.gz http://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz` (see header of `scripts/build_exon_index.R`) |
| `data/reference_exon_index.rds` | `Rscript scripts/build_exon_index.R` (needs the GTF above + `rtracklayer`; takes a few minutes) |
| `data/UP000005640_9606_reviewed_canonical.fasta`, `data/reference_mass_index.rds` | `Rscript scripts/build_reference_proteome.R` (downloads UniProt reviewed canonical proteome) |
| `data/reference_mz_index_{denatured,native}.rds` | `Rscript scripts/build_mz_collision_index.R` |
| `data/reference_digest_pool_<protease>.rds` | Built lazily on first middle-down "Run analysis" per protease (~20-40s, cached after) — no script needed, just use the app once per protease you need |
| `reference/genome/GRCh38.fa` + `.fai` | Download Ensembl's primary-assembly FASTA, gunzip, then `samtools faidx reference/genome/GRCh38.fa`:<br>`curl -o reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz http://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`<br>`gunzip -k reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` (rename result to `GRCh38.fa`)<br>`samtools faidx reference/genome/GRCh38.fa` |
| `reference/genome/GRCh38.mmi` | `minimap2 -d reference/genome/GRCh38.mmi reference/genome/GRCh38.fa` (large — ~12GB, takes a while) |
| `data/ensembl_exon_cache/`, `data/ensembl_protein_cache/` | Populated automatically on first use (live Ensembl REST calls, disk-cached per transcript — safe to delete and let it rebuild) |

All of the above use **relative paths from the project root** — nothing is
hardcoded to a machine-specific absolute path, so the project folder can be
moved/cloned anywhere without code changes.

## Directory structure

- `ui.R` / `server.R` / `global.R` — the Shiny app itself. `global.R` sources
  every `R/*.R` file and loads the precomputed indexes at startup.
- `R/proteoform_schema.R` — the shared `Proteoform` object every input module
  (gene/isoform selection, FASTA/ORF, rMATS) must resolve into before any
  downstream code (mass calc, resolvability scoring, fragment ladder) runs.
- `R/gene_isoform_data.R`, `R/ensembl_protein_fetch.R` — Option 1 (gene →
  isoform → proteoform), live Ensembl REST fetches with disk caching.
- `R/fasta_pipeline.R` — Option 2 (pasted FASTA → TransDecoder ORF →
  minimap2 spliced alignment → matched known transcripts).
- `R/rmats_adapter.R` — Option 3 (rMATS SE/MXE event → matched real
  transcripts + constructed synthetic isoforms). The most recently active and
  most subtle module — see `docs/DEVLOG.md` before modifying.
- `R/mass_calculation.R`, `R/resolving_power.R`, `R/charge_envelope.R`,
  `R/isotope_envelope.R`, `R/ms1_scoring.R` — the scoring engine (Phase 1-2).
- `R/fragment_ladder.R`, `R/fragmentation_propensity.R`,
  `R/fragment_collision.R` — MS2 fragment-ladder generation and cross-protein
  collision checks (Phase 3).
- `R/digestion.R` — middle-down protease-digestion simulation.
- `R/viz_json.R`, `R/svg_render.R`, `R/exon_axis_alignment.R`,
  `www/ptracker_viz.js` — visualization payload assembly (R side) and
  rendering (client-side JS; SVGs are built as HTML strings and injected via
  `session$sendCustomMessage()`, not React/htmlwidgets).
- `tests/testthat/` — run via `Rscript tests/testthat.R`.

## A few things that will save you time

- **Coordinates in rMATS files are 0-based-start** (`*ES` columns need `+1`);
  rMATS' own "upstream"/"downstream" column names are by **genomic**
  coordinate, not transcription direction — they flip relative to 5'/3' on a
  minus-strand gene. See `genomic_flanks_in_transcript_order()`'s doc comment.
- **`reference_exon_index` is CDS-only** (built from GTF `CDS` rows) — it has
  no UTR exons. Anywhere you need a transcript's *full* exon structure
  (UTRs included), use `fetch_transcript_exons()` (live Ensembl REST,
  disk-cached), not `reference_exon_index`.
- **Any `renderUI()` block whose own inputs it reads back can create a reset
  loop**: if widget A's value feeds a reactive that widget A's own `renderUI()`
  depends on, every change to A triggers a full re-render that recreates A —
  and if the recreated widget's `selected=`/`value=` isn't read back via
  `isolate(input$A)` first, the user's own pick gets silently discarded. This
  bit the rMATS backbone-transcript dropdown once already (see DEVLOG).
- Outputs inside a `conditionalPanel` that can become visible via a
  **server-side** switch (not the user's own click) need
  `outputOptions(output, "x", suspendWhenHidden = FALSE)` — otherwise Shiny
  never actually computes them on a fresh page load. Several outputs already
  have this; if you add a new one inside `force_gene_view`'s panel, add it
  there too.

## Testing

```r
Rscript tests/testthat.R
```
