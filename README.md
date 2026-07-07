# ProteoformTracker

Planning tool for isoform/proteoform-level detectability in top-down proteomics (TDP). Companion to [IsoPepTracker](https://github.com/HuangLabAtUAB/IsoPepTracker/) (bottom-up). R Shiny.

See `ptracker_vibe.md` (design spec) for full scope, scoring model, and roadmap.

## Status

**Phase 1 — core data model and mass engine** (in progress)
- [x] `Proteoform` schema (`R/proteoform_schema.R`)
- [x] Mass calculation with N/C-terminal processing + optional PTMs (`R/mass_calculation.R`)
- [x] Offline precomputed reference-proteome mass index (`R/reference_proteome_index.R`)

Phases 2-4 (resolvability scoring, fragmentation/MS2, validation) not yet started.

## Setup

Requires R >= 4.0 and Python (a project-local virtualenv is created automatically; the wrapped
`pyteomics` library is pinned below 5.0 -- see `python/requirements.txt` for why).

```r
# 1. Restore R package versions
renv::restore()

# 2. Create the Python venv + install pinned pyteomics
Rscript scripts/setup_python_env.R
```

Run the test suite:

```sh
Rscript tests/testthat.R
```

Run the app:

```r
shiny::runApp()
```

## Reference implementation

`reference/IsoPepTracker` (gitignored, local-only clone) holds the sibling project this tool
reuses code/conventions from -- GTF/GFF parsing, minimap2 wrapper, TransDecoder wrapper,
exon-track visualization, peptide-uniqueness core logic. Not part of this repo's history.
