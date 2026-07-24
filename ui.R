fluidPage(
  tags$head(
    tags$style(HTML("
    .pt-logo { font-size: 30px; font-weight: 700; margin-bottom: 2px; }
    .pt-tagline { color: #666; margin-bottom: 18px; }
    .pt-card { border: 1px solid #ddd; border-radius: 8px; padding: 14px 16px; margin-bottom: 10px; }
    .pt-card.active { border-color: #2a78d6; background: #f4f8fd; }
    .pt-card h4 { margin-top: 0; }
    .pt-settings { background: #fafafa; border: 1px solid #eee; border-radius: 8px; padding: 12px 16px; margin-bottom: 18px; }
    .pt-isorow { display: flex; align-items: center; gap: 10px; padding: 4px 0; border-bottom: 1px solid #eee; flex-wrap: wrap; }
    .pt-isorow .pt-id { font-family: monospace; min-width: 150px; }
    .pt-isorow .pt-meta { color: #777; font-size: 12px; min-width: 140px; }
    .pt-warn { color: #b03a2e; font-size: 12px; }
    .pt-note { color: #666; font-size: 12px; }
    svg.pt-viz { width: 100%; height: auto; border: 1px solid #eee; border-radius: 6px; cursor: default; }
    .pt-info-box { min-height: 20px; font-size: 12.5px; color: #333; padding: 6px 0; }
    .pt-zoom-controls { display: flex; align-items: center; gap: 6px; padding: 4px 0; }
    .pt-zoom-controls button { font-size: 13px; font-weight: 600; width: 26px; height: 26px; line-height: 1; border-radius: 5px; border: 1px solid #ccc; background: #f5f5f5; cursor: pointer; }
    .pt-zoom-controls button:hover { background: #e8e8e8; }
    .pt-zoom-controls .pt-zoom-reset { width: auto; padding: 0 8px; font-weight: 400; font-size: 11.5px; }
    .pt-zoom-controls .pt-zoom-range { font-size: 11px; color: #888; margin-left: 4px; }
    #pt-hover-tooltip { position: fixed; z-index: 10000; pointer-events: none; background: #222; color: #fff; font-size: 12px; line-height: 1.4; padding: 6px 9px; border-radius: 6px; box-shadow: 0 2px 8px rgba(0,0,0,.25); max-width: 260px; display: none; }
    .pt-mode-group { background: #fafafa; border: 1px solid #eee; border-radius: 8px; padding: 12px 16px; margin-bottom: 18px; }
    .pt-mode-group .radio-inline { font-size: 16px; font-weight: 600; color: #333; margin-right: 24px; }
    .pt-mode-group input[type='radio'] { transform: scale(1.25); margin-right: 8px; vertical-align: middle; }
    .pt-pep-table { width: 100%; border-collapse: collapse; font-size: 12.5px; }
    .pt-pep-table th, .pt-pep-table td { padding: 4px 8px; text-align: left; border-bottom: 1px solid #eee; white-space: nowrap; }
    .pt-pep-table th { color: #555; font-weight: 600; background: #fafafa; position: sticky; top: 0; }
    .pt-pep-table .form-group { margin-bottom: 0; }
    /* This app's Bootstrap (3.x) has no .modal-xl rule -- modalDialog(size='xl')
       silently falls back to the ~600px default width. Force a wide modal
       directly rather than depending on a size class the theme doesn't define. */
    .modal-dialog.modal-xl { width: 95vw; max-width: 1400px; }
    #btn_run_analysis:disabled { opacity: 0.5; cursor: not-allowed; }
    .pt-mode-group label.pt-disabled-label { opacity: 0.45; cursor: not-allowed; }
    #btn_load_gene:disabled, #btn_run_fasta:disabled, #btn_run_rmats:disabled,
    #btn_set_ms_strategy:disabled, #btn_set_ms_resolution:disabled { opacity: 0.5; cursor: not-allowed; }
  ")),
    tags$script(src = paste0("ptracker_viz.js?v=", as.integer(file.mtime("www/ptracker_viz.js"))))
  ),

  div(class = "pt-logo", "ProteoformTracker"),
  div(class = "pt-tagline", "Planning tool for isoform/proteoform-level detectability in top-down and middle-down proteomics"),

  div(class = "pt-settings",
    fluidRow(
      column(9, h4("MS strategy")),
      column(3, div(style = "text-align:right;padding-top:8px;",
        actionButton("btn_set_ms_strategy", "Set", class = "btn-default")))
    ),
    p(class = "pt-note", "Top-down analyzes the intact proteoform. Middle-down simulates a limited (partial) protease digestion first, then runs the same MS1/MS2 analysis on the resulting large peptides instead of the intact protein."),
    radioButtons("ms_strategy", NULL, inline = TRUE,
      choices = c("Top-down" = "topdown", "Middle-down" = "middledown"),
      selected = "topdown"
    ),
    conditionalPanel("input.ms_strategy == 'topdown'",
      fluidRow(
        column(3, numericInput("topdown_mass_min_kda", "Min protein mass (kDa)", value = 10, min = 0, step = 1)),
        column(3, numericInput("topdown_mass_max_kda", "Max protein mass (kDa)", value = 220, min = 1, step = 1))
      ),
      p(class = "pt-note", "Default 10-220 kDa covers the 2.5th-97.5th percentile (~95%) of the reviewed human proteome's intact monoisotopic mass, excluding both small fragments and the long tail of very large proteins (titin, dystrophin, etc.) that are impractical top-down MS targets. Isoforms outside this range are hidden from selection below, and the confounding-protein search pool is limited to it too -- adjust freely for your instrument's real usable mass range.")
    ),
    conditionalPanel("input.ms_strategy == 'middledown'",
      fluidRow(
        column(3, selectInput("md_protease", "Protease", choices = c("OmpT", "Lys-C", "Lys-N", "Glu-C", "Asp-N"))),
        column(3, numericInput("md_mass_min_kda", "Min peptide mass (kDa)", value = 3, min = 0.5, step = 0.5)),
        column(3, numericInput("md_mass_max_kda", "Max peptide mass (kDa)", value = 10, min = 1, step = 0.5))
      ),
      p(class = "pt-note", "Lys-C, Lys-N, and Glu-C are also used in bottom-up proteomics; under a limited (short-time) middle-down digestion they leave missed-cleavage sites, which this mass window is simulating rather than modeling digestion kinetics directly.")
    )
  ),

  div(class = "pt-settings",
    fluidRow(
      column(9, h4("MS resolution parameters")),
      column(3, div(style = "text-align:right;padding-top:8px;",
        actionButton("btn_set_ms_resolution", "Set", class = "btn-default")))
    ),
    p(class = "pt-note", "These feed every downstream scoring step (resolving power, envelope crowding, confounder search) regardless of which input path below you use."),
    fluidRow(
      column(3, numericInput("ms_r_ref", "Resolving power R (at reference m/z)", value = 120000, min = 1000, step = 10000)),
      column(3, numericInput("ms_mz_ref", "Reference m/z", value = 200, min = 50, step = 10)),
      column(3, numericInput("ms_safety_margin", "Safety margin (x FWHM)", value = 1.75, min = 1, max = 5, step = 0.05)),
      column(3, selectInput("ms_mode", "Ionization mode", choices = c("Denatured" = "denatured", "Native" = "native")))
    )
  ),

  div(class = "pt-mode-group",
    fluidRow(
      column(9, h4("Input selection")),
      column(3, div(style = "text-align:right;padding-top:8px;",
        actionButton("btn_reset_all", "Reset", class = "btn-default")))
    ),
    radioButtons("input_mode", NULL, inline = TRUE,
      choices = c(
        "1. Gene -> isoform -> proteoform" = "gene",
        "2. FASTA sequence" = "fasta",
        "3. rMATS alternative-splicing results" = "rmats"
      ),
      selected = "gene"
    ),
    # Hidden flag, distinct from input_mode itself: "Add to comparison" (in
    # module 2) needs to show the same merged proteoform/MS1/MS2 panel
    # module 1 uses (same analysis code either way), but flipping the
    # visible radio itself to "1. Gene..." reads as if the user had left
    # module 2, when the data underneath is still FASTA-derived. Setting
    # this instead (server.R) keeps the radio showing "2. FASTA sequence"
    # while still revealing the shared results panel below.
    tags$div(style = "display:none;", checkboxInput("force_gene_view", NULL, value = FALSE))
  ),
  hr(),

  # ---------------- Option 1: gene -> isoform -> proteoform ----------------
  conditionalPanel("input.input_mode == 'gene' || input.force_gene_view",
    div(class = "pt-card active",
      h4("Gene / isoform / proteoform selection"),
      fluidRow(
        column(4, textInput("gene_symbol", "Gene symbol", value = "", placeholder = "e.g. CD44")),
        column(2, br(), actionButton("btn_load_gene", "Load isoforms", class = "btn-primary"))
      ),
      textOutput("gene_status"),
      hr(),
      h5("Isoform catalog (real Ensembl transcripts for this gene)"),
      uiOutput("isoform_catalog_ui"),
      uiOutput("ptm_warnings_ui"),
      hr(),
      h5("Result proteoform table"),
      p(class = "pt-note", "Every checked isoform (plus its parsed PTM combinations) is included in the comparison below. Pick one as the confounder-search target:"),
      tableOutput("proteoform_table"),
      conditionalPanel("input.ms_strategy == 'middledown'",
        hr(),
        h5("Middle-down peptide candidates"),
        fluidRow(
          column(8, textOutput("digestion_summary_text")),
          column(4, actionButton("btn_open_peptide_picker", "Select peptides...", class = "btn-primary"))
        )
      ),
      selectInput("pf_target_select", "Confounder-search target", choices = character(0), width = "420px"),
      actionButton("btn_run_analysis", "Run analysis", class = "btn-success"),
      textOutput("stale_notice"),
      uiOutput("viz_script"),
      hr(),
      h5("1. Relevant-proteoform comparison"),
      p(class = "pt-note", "MS1 charge-envelope overlay, then the MS2 fragment ladder aligned on a shared exon axis. Scroll/pinch or use the +/- buttons to zoom, drag to pan. Hover a b/y tick once zoomed in enough to inspect it."),
      tags$svg(id = "s1-ms1", class = "pt-viz", viewBox = "0 0 640 130"),
      div(id = "s1-legend"),
      div(id = "s1-filters"),
      div(id = "s1-zoom"),
      tags$svg(id = "s1-ladder", class = "pt-viz", viewBox = "0 0 640 40"),
      div(id = "s1-info", class = "pt-info-box"),
      hr(),
      h5("2. Confounding-protein search (single target)"),
      textOutput("confounder_status"),
      tags$svg(id = "s2-ms1", class = "pt-viz", viewBox = "0 0 640 130"),
      div(id = "s2-legend"),
      div(id = "s2-filters"),
      div(id = "s2-zoom"),
      tags$svg(id = "s2-ladder", class = "pt-viz", viewBox = "0 0 640 40"),
      div(id = "s2-info", class = "pt-info-box")
    )
  ),

  # ---------------- Option 2: FASTA upload ----------------
  conditionalPanel("input.input_mode == 'fasta' && !input.force_gene_view",
    div(class = "pt-card",
      h4("FASTA sequence input"),
      p(class = "pt-note", "Provide a single spliced mRNA/cDNA (or CDS) nucleotide sequence -- not raw genomic DNA with introns. Two independent steps run: minimap2 spliced-aligns it against GRCh38 to identify which known gene/isoforms it belongs to (this does not depend on translation at all), and TransDecoder separately finds candidate open reading frames. Pick a known isoform to compare its exon structure against, pick which ORF candidate to use as the translation, add PTMs, then send it into the same comparison view as Option 1."),
      tabsetPanel(
        tabPanel("Paste sequence",
          textAreaInput("fasta_text", "Paste FASTA", rows = 6, placeholder = ">my_transcript\nACGTACGT...")
        ),
        tabPanel("Upload file",
          fileInput("fasta_file", "Choose FASTA file", accept = c(".fa", ".fasta", ".fas", ".txt"))
        )
      ),
      actionButton("btn_run_fasta", "Run sequence analysis", class = "btn-primary"),
      br(), br(),
      textOutput("fasta_status"),
      uiOutput("fasta_results_ui")
    )
  ),

  # ---------------- Option 3: rMATS upload ----------------
  conditionalPanel("input.input_mode == 'rmats' && !input.force_gene_view",
    div(class = "pt-card",
      h4("rMATS alternative-splicing results"),
      p(class = "pt-note", "rMATS only reports the differential exon(s) and their immediate flanking exons, not the rest of the transcript -- so a full-length proteoform can't be computed from the event alone. Instead, ProteoformTracker looks up which already-annotated transcripts of the gene (in the same precomputed exon index Option 1 uses) structurally match each arm of the event (exon-inclusion vs. exon-skipping for SE; 1st-exon vs. 2nd-exon for MXE), so you get real, full-length proteoforms rather than just the local differential region. If no annotated transcript matches an arm (some events reflect a splicing pattern no single annotated transcript uses), that arm shows no candidates -- constructing a synthetic transcript for that case isn't implemented yet. Only SE (skipped-exon) and MXE (mutually-exclusive-exons) events are supported so far; A3SS/A5SS/RI are a planned follow-up."),
      fluidRow(
        column(4, selectInput("rmats_event_type", "Event type", choices = c("SE", "MXE"))),
        column(6, fileInput("rmats_file", "Choose rMATS SE/MXE results file", accept = c(".txt", ".JC.txt")))
      ),
      textOutput("rmats_parse_status"),
      selectInput("rmats_event_select", "Event to analyze", choices = character(0), width = "560px"),
      actionButton("btn_run_rmats", "Find matching transcripts", class = "btn-primary"),
      br(), br(),
      textOutput("rmats_status"),
      uiOutput("rmats_results_ui"),
      # Static (not dynamically-generated) SVG/zoom/legend containers, same
      # as Option 1's s1-ladder/s1-ms1 -- unlike Option 2's exon-align
      # preview (whose containers live inside fasta_results_ui and are only
      # ever targeted by a message sent from an observer that necessarily
      # fires AFTER that UI exists, since it depends on widgets
      # fasta_results_ui itself creates), Option 3's preview message is sent
      # from the SAME observer that populates rmats_results_ui, racing that
      # renderUI's own client push -- keeping these elements always-present
      # in the DOM from page load avoids that race entirely (confirmed
      # necessary: with these inside rmats_results_ui instead, the custom
      # message reliably arrived before the SVG element existed and
      # silently did nothing).
      div(id = "rmats-exon-zoom"),
      tags$svg(id = "rmats-exon-align", class = "pt-viz", viewBox = "0 0 640 160"),
      div(id = "rmats-exon-align-legend")
    )
  )
)
