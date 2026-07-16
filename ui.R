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
  ")),
    tags$script(src = paste0("ptracker_viz.js?v=", as.integer(file.mtime("www/ptracker_viz.js"))))
  ),

  div(class = "pt-logo", "ProteoformTracker"),
  div(class = "pt-tagline", "Planning tool for isoform/proteoform-level detectability in top-down and middle-down proteomics"),

  div(class = "pt-settings",
    h4("MS strategy"),
    p(class = "pt-note", "Top-down analyzes the intact proteoform. Middle-down simulates a limited (partial) protease digestion first, then runs the same MS1/MS2 analysis on the resulting large peptides instead of the intact protein."),
    radioButtons("ms_strategy", NULL, inline = TRUE,
      choices = c("Top-down" = "topdown", "Middle-down" = "middledown"),
      selected = "topdown"
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
    h4("MS resolution parameters"),
    p(class = "pt-note", "These feed every downstream scoring step (resolving power, envelope crowding, confounder search) regardless of which input path below you use."),
    fluidRow(
      column(3, numericInput("ms_r_ref", "Resolving power R (at reference m/z)", value = 120000, min = 1000, step = 10000)),
      column(3, numericInput("ms_mz_ref", "Reference m/z", value = 200, min = 50, step = 10)),
      column(3, numericInput("ms_safety_margin", "Safety margin (x FWHM)", value = 1.75, min = 1, max = 5, step = 0.05)),
      column(3, selectInput("ms_mode", "Ionization mode", choices = c("Denatured" = "denatured", "Native" = "native")))
    )
  ),

  div(class = "pt-mode-group",
    h4("Input selection"),
    radioButtons("input_mode", NULL, inline = TRUE,
      choices = c(
        "1. Gene -> isoform -> proteoform" = "gene",
        "2. FASTA sequence" = "fasta",
        "3. rMATS alternative-splicing results" = "rmats"
      ),
      selected = "gene"
    )
  ),
  hr(),

  # ---------------- Option 1: gene -> isoform -> proteoform ----------------
  conditionalPanel("input.input_mode == 'gene'",
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
        p(class = "pt-note", "Every in-silico digest fragment (any number of missed cleavages) of the checked proteoforms above that falls in the mass window. Sorted by likely feasibility first (fewer missed cleavages, tighter MS1 peak), then by how many PTM sites it covers. Pick which candidate(s) to treat as \"proteins\" for the MS1/MS2 analysis below -- the top-ranked candidate per parent is pre-checked as a starting point."),
        uiOutput("digestion_candidates_ui"),
        textOutput("digestion_coverage_text")
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
  conditionalPanel("input.input_mode == 'fasta'",
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
      actionButton("btn_run_fasta", "Align (minimap2) + translate (TransDecoder)", class = "btn-primary"),
      br(), br(),
      textOutput("fasta_status"),
      uiOutput("fasta_results_ui")
    )
  ),

  # ---------------- Option 3: rMATS upload ----------------
  conditionalPanel("input.input_mode == 'rmats'",
    div(class = "pt-card",
      h4("rMATS alternative-splicing results"),
      p(class = "pt-note", "Upload one rMATS output file (SE/A3SS/A5SS/MXE/RI) to convert its AS events into updated exon structures for comparison."),
      selectInput("rmats_event_type", "Event type", choices = c("SE", "A3SS", "A5SS", "MXE", "RI")),
      fileInput("rmats_file", "Choose rMATS file", accept = c(".txt", ".JC.txt")),
      actionButton("btn_run_rmats", "Process rMATS file", class = "btn-primary"),
      br(), br(),
      textOutput("rmats_status")
    )
  )
)
