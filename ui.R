fluidPage(
  tags$head(
    tags$style(HTML("
    @font-face {
      font-family: 'Plex Mono';
      font-weight: 500;
      font-style: normal;
      src: url('fonts/plex-mono-500.woff2') format('woff2');
      font-display: swap;
    }
    @font-face {
      font-family: 'Plex Mono';
      font-weight: 700;
      font-style: normal;
      src: url('fonts/plex-mono-700.woff2') format('woff2');
      font-display: swap;
    }
    .pt-logo { font-family: 'Plex Mono', monospace; font-size: 30px; font-weight: 700; margin-bottom: 2px; letter-spacing: -0.01em; }
    .pt-tagline { font-family: 'Plex Mono', monospace; font-weight: 500; color: #666; margin-bottom: 18px; }
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
    /* Transcript-id sequence popover (hover ENST id in an isoform table to
       look up residue positions for PTM specs) -- deliberately a SEPARATE
       element from #pt-hover-tooltip above: that one is pointer-events:none
       (fine for a quick glance) and only 260px wide, both wrong for
       something meant to be read carefully and select/copy from. This one
       accepts pointer events (so moving the mouse INTO it keeps it open
       instead of it vanishing the instant you leave the trigger) and is
       wide enough for a monospace sequence block. */
    .pt-id.pt-seq-hover { cursor: help; border-bottom: 1px dotted #888; }
    #pt-seq-popover { position: fixed; z-index: 10001; background: #fff; color: #222; border: 1px solid #ccc; border-radius: 6px; box-shadow: 0 4px 16px rgba(0,0,0,.2); padding: 8px 10px; max-width: 640px; max-height: 420px; overflow: auto; display: none; }
    #pt-seq-popover .pt-seq-popover-header { font-weight: 600; font-size: 12.5px; margin-bottom: 4px; }
    #pt-seq-popover pre { margin: 0; font-family: monospace; font-size: 11px; line-height: 1.5; white-space: pre; }
    /* Collapsible isoform-catalog header (long isoform lists otherwise
       force scrolling all the way past a list the user is already done
       picking from, just to reach the results below). Triangle rotates
       -90deg when collapsed; the summary span shows a running selected-
       count so collapsing doesn't hide whether a pick actually stuck. */
    .pt-collapsible-header { cursor: pointer; display: flex; align-items: center; user-select: none; }
    .pt-collapsible-header:hover .pt-collapse-triangle { color: #2a78d6; }
    .pt-collapse-triangle { display: inline-block; font-size: 11px; color: #666; transition: transform 0.15s ease; }
    .pt-collapse-triangle.pt-collapsed { transform: rotate(-90deg); }
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
    /* MS1 aggregate stat tiles (clean/overlapping/total charge-state peaks),
       shown beside the MS1 chart itself (.pt-ms1-with-stats). Per-proteoform
       MS2 tallies are no longer a separate HTML panel at all -- they're
       drawn inline inside the ladder SVG next to each row's own title
       (rowBadgesSvg() in ptracker_viz.js), recomputed live against whichever
       fragment filter is currently active. */
    .pt-ms1-with-stats { display: flex; gap: 16px; align-items: flex-start; flex-wrap: wrap; }
    .pt-ms1-main { flex: 1 1 420px; min-width: 0; }
    .pt-stats-strip { flex: 0 0 190px; display: flex; flex-direction: column; gap: 8px; }
    .pt-stat-group-title { font-family: 'Plex Mono', monospace; font-size: 11px; font-weight: 700; color: #888; text-transform: uppercase; letter-spacing: 0.04em; }
    .pt-stat-tile { border: 1px solid #e2e2e2; border-radius: 8px; padding: 6px 10px; display: flex; align-items: center; gap: 8px; background: #fff; }
    .pt-stat-icon { font-size: 14px; width: 16px; text-align: center; flex-shrink: 0; display: inline-block; }
    .pt-stat-value { font-family: 'Plex Mono', monospace; font-weight: 700; font-size: 15px; line-height: 1.1; color: #222; }
    .pt-stat-label { font-size: 10.5px; color: #777; line-height: 1.2; }
    /* Visual separation between the MS1 overlay and the MS2 fragment-ladder
       section below it, so the two panels read as distinct rather than one
       continuous stack. */
    .pt-viz-section-heading { font-family: 'Plex Mono', monospace; font-size: 12.5px; font-weight: 700; color: #444; text-transform: uppercase; letter-spacing: 0.03em; margin: 18px 0 6px 0; padding-top: 14px; border-top: 1px solid #eee; }
    .pt-viz-section-heading:first-child { margin-top: 4px; padding-top: 0; border-top: none; }
    .pt-mode-group label.pt-disabled-label { opacity: 0.45; cursor: not-allowed; }
    #btn_load_gene:disabled, #btn_run_fasta:disabled, #btn_run_rmats:disabled,
    #btn_set_ms_strategy:disabled, #btn_set_ms_resolution:disabled { opacity: 0.5; cursor: not-allowed; }
  ")),
    tags$script(src = paste0("ptracker_viz.js?v=", as.integer(file.mtime("www/ptracker_viz.js"))))
  ),

  div(style = "display:flex; align-items:center; gap:14px; margin-bottom:18px;",
    HTML(r"(<svg width="180" height="100" viewBox="36 118 409 228" role="img" aria-label="ProteoformTracker icon">
<rect x="56" y="181" width="36" height="36" rx="4" fill="#1f9c86"/>
<rect x="112" y="181" width="40" height="36" rx="4" fill="#1f9c86"/>
<rect x="172" y="181" width="36" height="36" rx="4" fill="#1f9c86"/>
<rect x="228" y="181" width="30" height="36" rx="4" fill="#1f9c86"/>
<line x1="74" y1="181" x2="74" y2="156" stroke="#ffd23f" stroke-width="4" stroke-linecap="round"/>
<circle cx="74" cy="149" r="11" fill="#ffd23f" stroke="#17332e" stroke-width="2.5"/>
<rect x="56" y="281" width="36" height="36" rx="4" fill="#ff6b4a"/>
<rect x="112" y="281" width="40" height="36" rx="4" fill="#ff6b4a"/>
<rect x="172" y="281" width="36" height="36" rx="4" fill="none" stroke="#b0b0b0" stroke-width="2" stroke-dasharray="4,4"/>
<rect x="228" y="281" width="30" height="36" rx="4" fill="#ff6b4a"/>
<path d="M300,273 C314,273 326,183 340,183 C354,183 366,273 380,273 Z" fill="#1f9c86" fill-opacity="0.8"/>
<path d="M345,273 C359,273 371,193 385,193 C399,193 411,273 425,273 Z" fill="#ff6b4a" fill-opacity="0.8"/>
<g stroke-linecap="round">
  <g stroke-width="3.5">
    <line x1="305" y1="288" x2="305" y2="298" stroke="#8a8a8a"/>
    <line x1="314" y1="288" x2="314" y2="303" stroke="#8a8a8a"/>
    <line x1="323" y1="288" x2="323" y2="300" stroke="#1f9c86"/>
    <line x1="332" y1="288" x2="332" y2="296" stroke="#8a8a8a"/>
    <line x1="350" y1="288" x2="350" y2="299" stroke="#8a8a8a"/>
    <line x1="359" y1="288" x2="359" y2="304" stroke="#8a8a8a"/>
    <line x1="368" y1="288" x2="368" y2="297" stroke="#1f9c86"/>
    <line x1="377" y1="288" x2="377" y2="301" stroke="#8a8a8a"/>
    <line x1="395" y1="288" x2="395" y2="298" stroke="#8a8a8a"/>
    <line x1="404" y1="288" x2="404" y2="302" stroke="#8a8a8a"/>
    <line x1="413" y1="288" x2="413" y2="299" stroke="#1f9c86"/>
  </g>
  <g stroke-width="3.5">
    <line x1="305" y1="310" x2="305" y2="320" stroke="#8a8a8a"/>
    <line x1="314" y1="310" x2="314" y2="325" stroke="#8a8a8a"/>
    <line x1="332" y1="310" x2="332" y2="318" stroke="#8a8a8a"/>
    <line x1="341" y1="310" x2="341" y2="324" stroke="#ff6b4a"/>
    <line x1="350" y1="310" x2="350" y2="321" stroke="#8a8a8a"/>
    <line x1="359" y1="310" x2="359" y2="326" stroke="#8a8a8a"/>
    <line x1="377" y1="310" x2="377" y2="323" stroke="#8a8a8a"/>
    <line x1="386" y1="310" x2="386" y2="322" stroke="#ff6b4a"/>
    <line x1="395" y1="310" x2="395" y2="320" stroke="#8a8a8a"/>
    <line x1="404" y1="310" x2="404" y2="324" stroke="#8a8a8a"/>
  </g>
</g>
</svg>)"),
    div(
      div(class = "pt-logo", style = "margin-bottom: 0;", "ProteoformTracker"),
      div(class = "pt-tagline", style = "margin-bottom: 0; font-size: 11.5px; letter-spacing: 0.05em; text-transform: uppercase;",
          "Detectability planning for top/middle-down proteomics"),
      div(class = "pt-note", style = "margin-top: 4px;",
          "Tutorial & documentation: ",
          tags$a(href = "https://huanglabatuab.github.io/proteoformtracker-docs/",
                 target = "_blank", rel = "noopener noreferrer",
                 "huanglabatuab.github.io/proteoformtracker-docs"))
    )
  ),

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
        column(4, textInput("gene_symbol", "Gene symbol", value = "", placeholder = "e.g. BCL2L1")),
        column(2, br(), actionButton("btn_load_gene", "Load isoforms", class = "btn-primary"))
      ),
      textOutput("gene_status"),
      hr(),
      tags$div(id = "isoform-catalog-toggle", class = "pt-collapsible-header",
        tags$span(id = "isoform-catalog-triangle", class = "pt-collapse-triangle", HTML("&#9660;")),
        h5(style = "display:inline; margin:0 0 0 6px;", "Isoform catalog (real Ensembl transcripts for this gene)"),
        tags$span(id = "isoform-catalog-summary", class = "pt-note", style = "margin-left:10px;")
      ),
      tags$div(id = "isoform-catalog-collapse-body",
        uiOutput("isoform_catalog_ui"),
        uiOutput("ptm_warnings_ui")
      ),
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
      tags$div(id = "s1-toggle", class = "pt-collapsible-header", `data-collapse-target` = "s1-collapse-body",
        tags$span(class = "pt-collapse-triangle pt-collapsed", HTML("&#9660;")),
        h5(style = "display:inline; margin:0 0 0 6px;", "1. Relevant-proteoform comparison"),
        tags$span(id = "s1-collapse-hint", class = "pt-note", style = "margin-left:10px;", "(run analysis to populate)")
      ),
      tags$div(id = "s1-collapse-body", style = "display:none;",
        div(class = "pt-viz-section-heading", "MS1 charge-envelope overlay"),
        p(class = "pt-note", "Scroll/pinch or use the +/- buttons to zoom, drag to pan -- individual isotope peaks are only visible once zoomed into a single charge state's own narrow m/z window; a smooth envelope is shown where the instrument wouldn't resolve them."),
        tags$div(class = "pt-ms1-with-stats",
          tags$div(class = "pt-ms1-main",
            tags$svg(id = "s1-ms1", class = "pt-viz", viewBox = "0 0 640 130"),
            div(id = "s1-ms1-zoom")
          ),
          tags$div(id = "s1-stats-strip", class = "pt-stats-strip")
        ),
        div(class = "pt-viz-section-heading", "MS2 fragment ladder"),
        div(class = "pt-mode-group", style = "margin-bottom:14px;",
          radioButtons("scoring_mode", "Fragmentation scoring mode", inline = TRUE,
            choices = c(
              "RF ranking (no length -- for long proteoforms)" = "rf",
              "Calibrated (length-aware)" = "glm"
            ),
            selected = "rf"
          )
        ),
        p(class = "pt-note", "Aligned on a shared exon axis. Scroll/pinch or use the +/- buttons to zoom, drag to pan. Hover a b/y tick once zoomed in enough to inspect it; the unique/partial/common counts next to each proteoform's title update live as you change the filter below."),
        downloadButton("btn_download_s1_peaks", "Download peak data (TSV)", class = "btn-default btn-sm", style = "margin-bottom:10px;"),
        div(id = "s1-legend"),
        div(id = "s1-filters"),
        div(id = "s1-zoom"),
        tags$svg(id = "s1-ladder", class = "pt-viz", viewBox = "0 0 640 40"),
        div(id = "s1-info", class = "pt-info-box"),
        p(class = "pt-note", "Click a b/y tick above (Elevated/High/Very high filter -- propensity score>1) to see that specific fragment ion's own isotope peaks below, overlaid across every checked proteoform that has a qualifying fragment at the same aligned position -- the same resolved-vs-overlapping question the MS1 chart above answers for the intact protein, answered here for one fragment ion at a time."),
        div(id = "s1-frag-label", class = "pt-note", style = "font-weight:600;"),
        tags$svg(id = "s1-frag-ms1", class = "pt-viz", viewBox = "0 0 640 130"),
        div(id = "s1-frag-ms1-zoom")
      ),
      hr(),
      tags$div(id = "s2-toggle", class = "pt-collapsible-header", `data-collapse-target` = "s2-collapse-body",
        tags$span(class = "pt-collapse-triangle pt-collapsed", HTML("&#9660;")),
        h5(style = "display:inline; margin:0 0 0 6px;", "2. Confounding-protein search (single target)"),
        tags$span(id = "s2-collapse-hint", class = "pt-note", style = "margin-left:10px;", "(run analysis to populate)")
      ),
      tags$div(id = "s2-collapse-body", style = "display:none;",
        uiOutput("confounder_status"),
        uiOutput("confounder_candidate_list_ui"),
        hr(),
        uiOutput("viz_script_s2"),
        tags$div(id = "s2-ms1-toggle", class = "pt-collapsible-header", `data-collapse-target` = "s2-ms1-collapse-body",
          tags$span(class = "pt-collapse-triangle pt-collapsed", HTML("&#9660;")),
          h5(style = "display:inline; margin:0 0 0 6px;", "MS1 charge-envelope overlay")
        ),
        tags$div(id = "s2-ms1-collapse-body", style = "display:none;",
          tags$div(class = "pt-ms1-with-stats",
            tags$div(class = "pt-ms1-main",
              tags$svg(id = "s2-ms1", class = "pt-viz", viewBox = "0 0 640 130"),
              div(id = "s2-ms1-zoom")
            ),
            tags$div(id = "s2-stats-strip", class = "pt-stats-strip")
          )
        ),
        div(id = "s2-legend"),
        tags$div(id = "s2-ms2-toggle", class = "pt-collapsible-header", `data-collapse-target` = "s2-ms2-collapse-body",
          tags$span(class = "pt-collapse-triangle pt-collapsed", HTML("&#9660;")),
          h5(style = "display:inline; margin:0 0 0 6px;", "MS2 fragment ladder")
        ),
        tags$div(id = "s2-ms2-collapse-body", style = "display:none;",
          downloadButton("btn_download_s2_peaks", "Download peak data (TSV)", class = "btn-default btn-sm", style = "margin-bottom:10px;"),
          div(id = "s2-filters"),
          div(id = "s2-zoom"),
          tags$svg(id = "s2-ladder", class = "pt-viz", viewBox = "0 0 640 40"),
          div(id = "s2-info", class = "pt-info-box"),
          p(class = "pt-note", "Click a b/y tick above (propensity score>1) to see that fragment ion's own isotope peaks below."),
          div(id = "s2-frag-label", class = "pt-note", style = "font-weight:600;"),
          tags$svg(id = "s2-frag-ms1", class = "pt-viz", viewBox = "0 0 640 130"),
          div(id = "s2-frag-ms1-zoom")
        )
      )
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
      p(class = "pt-note", "rMATS only reports the differential exon(s) and their immediate flanking exons, not the rest of the transcript -- so a full-length proteoform can't be computed from the event alone. Instead, ProteoformTracker looks up which already-annotated transcripts of the gene (in the same precomputed exon index Option 1 uses) structurally match each arm of the event (exon-inclusion vs. exon-skipping for SE; 1st-exon vs. 2nd-exon for MXE; intron-retained vs. -spliced for RI; long- vs. short-exon form for A5SS/A3SS), so you get real, full-length proteoforms rather than just the local differential region. Every arm also gets a \"constructed\" synthetic isoform (a user-pickable backbone transcript with the local region replaced by rMATS' own reported exons), so an arm with no real annotated match still has something usable."),
      fluidRow(
        column(4, selectInput("rmats_event_type", "Event type", choices = c("SE", "MXE", "RI", "A5SS", "A3SS"))),
        column(6, fileInput("rmats_file", "Choose rMATS results file", accept = c(".txt", ".JC.txt")))
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
