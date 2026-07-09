# ProteoformTracker

Planning tool for isoform/proteoform-level detectability in top-down proteomics (TDP). Companion to [IsoPepTracker](https://github.com/HuangLabAtUAB/IsoPepTracker/) (bottom-up). R Shiny.

See `ptracker_vibe.md` (design spec) for full scope, scoring model, and roadmap.

## Status

**Phase 1 — core data model and mass engine** (done)
- [x] `Proteoform` schema (`R/proteoform_schema.R`), including N-terminal Met excision and
      user-supplied signal peptide/propeptide cleavage sites
- [x] Mass calculation with N/C-terminal processing + optional PTMs (`R/mass_calculation.R`)
- [x] Curated Unimod lookup table for common PTMs (`R/unimod_table.R`, `unimod_ptm("Phospho", site)`) --
      not the full Unimod database; falls back to `ptm()` with a manual mass delta for anything else
- [x] Offline precomputed reference-proteome mass index (`R/reference_proteome_index.R`), built from
      the human reviewed canonical proteome (UniProt `UP000005640`, 20,391 sequences after filtering
      non-standard residues) via `scripts/build_reference_proteome.R`

**Phase 2 — resolvability scoring and MS1 visualization** (scoring engine done; visualization not yet wired into the app)
- [x] Orbitrap R(m/z) resolving-power model + mass-domain ΔM_FWHM (`R/resolving_power.R`)
- [x] Charge-state envelope prediction (basic-residue-count ceiling for small denatured proteoforms,
      Rayleigh-limit scaling otherwise) + relative S/N penalty model (`R/charge_envelope.R`)
- [x] Averagine isotope-envelope width model (`R/isotope_envelope.R`) -- closed-form binomial-moments
      approximation, not full combinatorial isotopologue enumeration (intractable at intact-protein scale)
- [x] MS1 resolvability verdict + envelope-crowding cross-proteoform m/z collision check
      (`R/ms1_scoring.R`)
- [ ] MS1 m/z envelope overlay visualization in the app (prototyped standalone, not yet wired into `ui.R`/`server.R`)

Phases 3-4 (fragmentation/MS2, validation) not yet started.

## Setup

Requires R >= 4.0 and Python (a project-local virtualenv is created automatically; the wrapped
`pyteomics` library is pinned below 5.0 -- see `python/requirements.txt` for why).

```r
# 1. Restore R package versions
renv::restore()

# 2. Create the Python venv + install pinned pyteomics
Rscript scripts/setup_python_env.R

# 3. Download the reference proteome and build the confounding-protein mass index
# (data/ is gitignored -- this regenerates it; takes ~1-2 min)
Rscript scripts/build_reference_proteome.R
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
