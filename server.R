function(input, output, session) {

  # ============================================================
  # Option 1: gene -> isoform -> proteoform
  # ============================================================

  catalog <- reactiveVal(NULL)
  bare_pfs <- reactiveVal(NULL)
  # Transcripts that should default to checked (beyond the "NOVEL_" prefix
  # rule), and a PTM-spec prefill, both set by Option 2's "Add to
  # comparison" step -- read at isoform_catalog_ui render time rather than
  # via updateCheckboxInput/updateTextInput, since those would race the
  # renderUI that first creates the widgets they'd target.
  precheck_extra <- reactiveVal(character(0))
  ptm_prefill <- reactiveVal(list())
  # Named list, representative transcript id -> every transcript id
  # (including itself) that shares its exact protein sequence -- see
  # dedupe_proteoforms_by_sequence(). catalog()$transcript_ids only ever
  # holds representatives; this is purely for the "+N more" display note.
  isoform_synonyms <- reactiveVal(list())
  # Which module last populated catalog()/bare_pfs(): "gene" | "fasta" |
  # "rmats" | NULL (nothing loaded / since Reset). Module 2 and 3 both
  # merge their results into the SAME catalog()/bare_pfs() Option 1 uses
  # (so the rest of the app -- proteoform table, digestion, Run analysis --
  # works unchanged regardless of where the data came from), but that means
  # "is Option 1's own Load-isoforms button safe to click" can't be told
  # apart from proteoform id naming alone (module 3's ids are real ENST
  # ids, same as module 1's, unlike module 2's "NOVEL_" prefix) -- this
  # flag is the general, source-agnostic way to know.
  catalog_source <- reactiveVal(NULL)
  # TRUE from the moment "Run analysis" or "Run sequence analysis" is
  # clicked until the user clicks Reset -- while TRUE, the 3 "Input
  # selection" radio choices are greyed out and unclickable, so the user
  # can't switch modules mid-analysis and get distracted by unrelated
  # module state instead of the results they just asked for.
  analysis_locked <- reactiveVal(FALSE)
  observe({
    session$sendCustomMessage("pt_set_selector_enabled",
      list(selector = 'input[name="input_mode"]', enabled = !analysis_locked()))
  })
  # TRUE only once "Run analysis" has actually completed for whatever is
  # CURRENTLY loaded; FALSE after Reset. Needed because `analysis` below is
  # an eventReactive -- it has no way to be told "return empty" and simply
  # keeps returning its last-computed value forever, frozen, until the next
  # click. Guarding viz_script/confounder_status on catalog() alone isn't
  # enough: Reset sets catalog(NULL) (correctly blanking the panel for that
  # instant), but the next gene loaded or FASTA "Add to comparison" makes
  # catalog() non-NULL again -- and without this flag, that alone was enough
  # to pass viz_script's guard and replay the STALE pre-reset analysis()
  # result straight into the DOM, before the user ever clicked "Run
  # analysis" again for the new data (reproduced directly: reset after a
  # TP53 top-down run, then loaded an RBMS1 FASTA sequence + "Add to
  # comparison" -- the ladder panel came back showing the old TP53 result).
  analysis_ready <- reactiveVal(FALSE)

  # "Set" buttons for the two global settings panels (MS strategy, MS
  # resolution parameters): clicking Set greys the button out to show the
  # user their choices there are locked in; changing anything in that same
  # panel un-sets it again, requiring a fresh explicit Set click. Load
  # isoforms / Run sequence analysis / Run analysis all additionally
  # require BOTH to be TRUE (see their own gating below) -- these are
  # deliberately session-wide, not per-analysis, so Reset does not touch
  # them (same reasoning already applied to the settings themselves).
  ms_strategy_set <- reactiveVal(FALSE)
  ms_resolution_set <- reactiveVal(FALSE)
  observeEvent(input$btn_set_ms_strategy, { ms_strategy_set(TRUE) })
  observeEvent(input$btn_set_ms_resolution, { ms_resolution_set(TRUE) })
  # Un-set on ANY change within that panel -- "mass filter" here covers both
  # the top-down intact-protein range and the middle-down protease/peptide
  # range, since which one is even relevant depends on ms_strategy itself.
  observeEvent(list(input$ms_strategy, input$topdown_mass_min_kda, input$topdown_mass_max_kda,
                     input$md_protease, input$md_mass_min_kda, input$md_mass_max_kda),
    { ms_strategy_set(FALSE) }, ignoreInit = TRUE)
  observeEvent(list(input$ms_r_ref, input$ms_mz_ref, input$ms_safety_margin, input$ms_mode),
    { ms_resolution_set(FALSE) }, ignoreInit = TRUE)
  observe({
    session$sendCustomMessage("pt_set_button_enabled", list(id = "btn_set_ms_strategy", enabled = !ms_strategy_set()))
    session$sendCustomMessage("pt_set_button_enabled", list(id = "btn_set_ms_resolution", enabled = !ms_resolution_set()))
  })

  # Format-only gene-symbol validation for gating "Load isoforms" -- real
  # gene lookup happens on click (build_gene_isoform_catalog()); this is
  # just "is there non-trivial, gene-symbol-shaped text here", cheap enough
  # to run on every keystroke with no network/lookup cost.
  is_valid_gene_symbol <- function(x) {
    x <- trimws(x %||% "")
    nzchar(x) && grepl("^[A-Za-z0-9_.-]+$", x)
  }

  # Format-only FASTA-content validation for gating "Run sequence analysis"
  # -- strips header line(s) and non-letter characters, then requires a
  # real run of nucleotide letters of a sane minimum length. Deliberately
  # not a full parser: just enough to distinguish "nothing/garbage pasted
  # or uploaded yet" from "there's a real sequence here".
  looks_like_valid_fasta <- function(text) {
    if (is.null(text) || !nzchar(trimws(text))) return(FALSE)
    lines <- strsplit(text, "\n")[[1]]
    seq_lines <- lines[!grepl("^>", lines)]
    seq <- toupper(gsub("[^A-Za-z]", "", paste(seq_lines, collapse = "")))
    nchar(seq) >= 30 && grepl("^[ACGTUN]+$", seq)
  }

  # Outputs that live inside the "gene" conditionalPanel get SUSPENDED by
  # default while that panel is hidden (e.g. while the user is on the
  # "fasta" tab) -- Shiny skips computing hidden outputs to save work. That
  # silently broke Option 2's "Add to comparison": it sets these render
  # functions and switches to the gene tab via updateRadioButtons() in the
  # SAME observer, but since the panel was still hidden at assignment time,
  # the new render never actually got evaluated/pushed to the client (found
  # by comparing Shiny.shinyapp.$values in the browser against what the
  # server had genuinely computed -- catalog()/bare_pfs() were correct
  # server-side the whole time; only the client-visible output was stuck on
  # its stale pre-hidden value). outputOptions() requires the output to
  # already exist, hence the placeholder assignments before it; the real
  # render functions defined later simply replace these.
  output$gene_status <- renderText("")
  output$isoform_catalog_ui <- renderUI(NULL)
  output$ptm_warnings_ui <- renderUI(NULL)
  output$proteoform_table <- renderTable(NULL)
  output$stale_notice <- renderText("")
  output$viz_script <- renderUI(NULL)
  output$confounder_status <- renderText("")
  # digestion_candidates_ui/digestion_coverage_text sit inside a SECOND,
  # nested conditionalPanel (ms_strategy == 'middledown' inside input_mode
  # == 'gene') -- same suspend-while-hidden gap as above: becoming visible
  # only via the user's own click (not a same-observer server-side tab
  # switch, unlike the original bug) still left these two permanently stuck
  # unrendered client-side (confirmed via Shiny.shinyapp.$values never
  # having the key at all), so they need the same placeholder + explicit
  # suspendWhenHidden = FALSE treatment. They now live inside a modalDialog
  # (opened on demand) rather than that conditionalPanel directly, which is
  # an even harder case -- their container doesn't exist in the DOM at all
  # until showModal() runs -- so keeping them un-suspended, so content is
  # already computed and ready the moment the modal opens, still applies.
  output$digestion_candidates_ui <- renderUI(NULL)
  output$digestion_coverage_text <- renderText("")
  output$digestion_summary_text <- renderText("")
  output$digestion_per_parent_ui <- renderUI(NULL)
  # fasta_status/fasta_results_ui (module 2) and rmats_status/
  # rmats_results_ui (module 3) hit this exact same gap: reproduced
  # directly for module 3 on a fresh page load (switch straight to the
  # rMATS tab, upload, click "Find matching transcripts" -- the observer
  # runs and assigns output$rmats_status/output$rmats_results_ui with no
  # error, but Shiny's client-side hidden-tracking for that conditionalPanel
  # never flips to "visible" on a page that never rendered any other tab
  # first, so the assignment is silently never pushed to the client).
  # Fixing module 2's pair too defensively, since nothing here is
  # module-3-specific -- module 2 just happened not to hit it in prior
  # ad-hoc testing, not because it's actually immune.
  output$fasta_status <- renderText("")
  output$fasta_results_ui <- renderUI(NULL)
  output$rmats_status <- renderText("")
  output$rmats_results_ui <- renderUI(NULL)
  for (nm in c("gene_status", "isoform_catalog_ui", "ptm_warnings_ui", "proteoform_table",
               "stale_notice", "viz_script", "confounder_status",
               "digestion_candidates_ui", "digestion_coverage_text", "digestion_summary_text",
               "digestion_per_parent_ui", "fasta_status", "fasta_results_ui",
               "rmats_status", "rmats_results_ui")) {
    outputOptions(output, nm, suspendWhenHidden = FALSE)
  }

  # Module 1's own gene_symbol/"Load isoforms" are disabled (and the box
  # cleared) whenever module 2 or 3 is in play -- either the user is
  # actively on the FASTA/rMATS tab, or they've already run "Add to
  # comparison" (which jumps the view back to the "gene" tab to show the
  # merged proteoform table, but bare_pfs() still holds the module-2/3-
  # derived entries it added, tracked via catalog_source()). Two real
  # problems this fixes: a stale gene symbol from an earlier, unrelated
  # module-1 session was left sitting in the box after switching modules,
  # reading as if it were still connected to the current data; and
  # clicking "Load isoforms" while a module-2/3 result is active would
  # silently overwrite catalog()/bare_pfs() entirely, destroying it with no
  # warning.
  observe({
    other_module_active <- input$input_mode %in% c("fasta", "rmats") ||
      identical(catalog_source(), "fasta") || identical(catalog_source(), "rmats")
    if (other_module_active) updateTextInput(session, "gene_symbol", value = "")
    load_gene_enabled <- !other_module_active && is_valid_gene_symbol(input$gene_symbol) &&
      ms_strategy_set() && ms_resolution_set()
    session$sendCustomMessage("pt_set_button_enabled", list(id = "btn_load_gene", enabled = load_gene_enabled))
    session$sendCustomMessage("pt_set_button_enabled", list(id = "gene_symbol", enabled = !other_module_active))
  })

  # force_gene_view (see ui.R) only makes sense as an override while the
  # radio itself still says "fasta"/"rmats" -- once the user explicitly
  # clicks a DIFFERENT tab themselves, drop the override so it doesn't
  # linger and show the wrong panel the next time they click back. ignoreInit:
  # this fires on the user's OWN radio clicks only, not the initial page
  # load (force_gene_view already starts FALSE).
  observeEvent(input$input_mode, {
    updateCheckboxInput(session, "force_gene_view", value = FALSE)
  }, ignoreInit = TRUE)

  observeEvent(input$btn_load_gene, {
    req(input$gene_symbol)
    gene <- toupper(trimws(input$gene_symbol))
    if (is.null(reference_exon_index)) {
      output$gene_status <- renderText("No exon index loaded -- run scripts/build_exon_index.R first.")
      return()
    }
    cat_result <- build_gene_isoform_catalog(gene, reference_exon_index)
    if (is.null(cat_result)) {
      output$gene_status <- renderText(sprintf('"%s" not found in the precomputed exon index.', gene))
      catalog(NULL)
      bare_pfs(NULL)
      isoform_synonyms(list())
      catalog_source(NULL)
      return()
    }

    n_transcripts <- length(cat_result$transcript_ids)
    pfs <- withProgress(message = paste("Fetching real protein sequences for", gene), value = 0, {
      n <- n_transcripts
      pfs <- list()
      fetch_one <- function(tid) {
        tryCatch({
          seq <- fetch_transcript_protein(tid)
          if (is.na(seq)) NULL else proteoform(id = tid, sequence = seq, provenance = "module1_isoform_selection")
        }, error = function(e) NULL)
      }
      for (i in seq_along(cat_result$transcript_ids)) {
        tid <- cat_result$transcript_ids[i]
        incProgress(1 / n, detail = tid)
        pfs[[tid]] <- fetch_one(tid)
      }
      # Ensembl's public REST API intermittently 500s/times out under a
      # burst of ~20-40 sequential requests (observed directly), even with
      # per-request retries -- a second pass just on whatever failed, after
      # a short pause, resolves most of those stragglers without slowing
      # down the common case where everything succeeds first try.
      failed_tids <- names(pfs)[vapply(pfs, is.null, logical(1))]
      if (length(failed_tids) > 0) {
        Sys.sleep(2)
        for (tid in failed_tids) pfs[[tid]] <- fetch_one(tid)
      }
      pfs
    })

    n_ok <- sum(!vapply(pfs, is.null, logical(1)))
    # Collapse transcripts that translate to an IDENTICAL protein sequence
    # down to one representative each (see dedupe_proteoforms_by_sequence()'s
    # doc comment for why this has to compare actual sequences rather than
    # Ensembl protein/ENSP accessions) -- shortens what's often a long,
    # heavily redundant transcript list down to the real, distinct set of
    # proteins, which is what actually matters for MS1/MS2 planning.
    ok_pfs <- pfs[!vapply(pfs, is.null, logical(1))]
    deduped <- dedupe_proteoforms_by_sequence(ok_pfs)
    bare_pfs(deduped$pf_list)
    isoform_synonyms(deduped$synonyms)

    rep_ids <- names(deduped$pf_list)
    cat_result$transcript_ids <- rep_ids
    cat_result$protein_lengths <- cat_result$protein_lengths[rep_ids]
    catalog(cat_result)
    catalog_source("gene")

    n_unique <- length(rep_ids)
    output$gene_status <- renderText(sprintf(
      "%s: %d transcripts in the precomputed exon index, %d with a real fetched protein sequence, %s.",
      gene, n_transcripts, n_ok,
      if (n_unique < n_ok) sprintf("collapsed to %d unique protein sequence(s)", n_unique) else sprintf("%d unique protein sequence(s)", n_unique)
    ))
  })

  output$isoform_catalog_ui <- renderUI({
    cat_result <- catalog()
    pfs <- bare_pfs()
    if (is.null(cat_result) || is.null(pfs)) return(tags$p(class = "pt-note", "Load a gene to see its real isoform catalog."))

    # Top-down mass filter (see "MS strategy" panel) only makes sense for
    # top-down: it's the intact protein that has to fall in the MS
    # instrument's usable mass range. Middle-down instead applies its own
    # mass window to the DIGESTED PEPTIDES later (digestion_result()), so a
    # large protein is still fine there even if its own intact mass would
    # be out of this range.
    topdown_filter <- identical(input$ms_strategy, "topdown")
    mass_min_da <- (input$topdown_mass_min_kda %||% 10) * 1000
    mass_max_da <- (input$topdown_mass_max_kda %||% 220) * 1000
    visible_ids <- if (!topdown_filter) cat_result$transcript_ids else Filter(function(tid) {
      pf <- pfs[[tid]]
      is.null(pf) || (proteoform_mass(pf)$mass >= mass_min_da && proteoform_mass(pf)$mass <= mass_max_da)
    }, cat_result$transcript_ids)
    n_hidden <- length(cat_result$transcript_ids) - length(visible_ids)

    rows <- lapply(visible_ids, function(tid) {
      pf <- pfs[[tid]]
      if (is.null(pf)) {
        return(tags$div(class = "pt-isorow",
          tags$input(type = "checkbox", disabled = NA),
          tags$span(class = "pt-id", tid),
          tags$span(class = "pt-meta", "no protein sequence available")
        ))
      }
      exon_nums <- transcript_exon_numbers(cat_result$exon_table, tid)
      mass <- proteoform_mass(pf)$mass
      prev_checked <- isolate(input[[paste0("iso_chk_", tid)]])
      default_checked <- if (!is.null(prev_checked)) prev_checked else (startsWith(tid, "NOVEL_") || tid %in% isolate(precheck_extra()))
      prev_ptm <- isolate(input[[paste0("iso_ptm_", tid)]])
      ptm_default <- if (!is.null(prev_ptm)) prev_ptm else (isolate(ptm_prefill())[[tid]] %||% "")
      # Transcripts sharing an identical protein sequence (different ENSP
      # per transcript regardless -- see dedupe_proteoforms_by_sequence())
      # are collapsed to this one representative row; note the others here
      # rather than silently dropping them from view.
      synonym_group <- isolate(isoform_synonyms())[[tid]]
      synonym_note <- if (!is.null(synonym_group) && length(synonym_group) > 1) {
        sprintf(" (+%d more transcript(s) with this same sequence: %s)",
                length(synonym_group) - 1, paste(setdiff(synonym_group, tid), collapse = ", "))
      } else ""
      tags$div(class = "pt-isorow",
        checkboxInput(paste0("iso_chk_", tid), NULL, value = default_checked),
        tags$span(class = "pt-id", tid),
        tags$span(class = "pt-meta", sprintf("%d aa, %.1f Da", nchar(pf$sequence), mass)),
        tags$span(class = "pt-meta", sprintf("exons %s", compress_exon_ranges(exon_nums))),
        textInput(paste0("iso_ptm_", tid), NULL, value = ptm_default,
                  placeholder = "e.g. 133_Thr_Phospho; 210_Pro_Oxidation,215_Ser_Sulfo", width = "420px"),
        if (nzchar(synonym_note)) tags$span(class = "pt-note", style = "font-size:11px;", synonym_note)
      )
    })
    tagList(
      if (topdown_filter && n_hidden > 0) tags$p(class = "pt-note",
        sprintf("%d isoform(s) outside the %g-%g kDa top-down mass filter are hidden from selection.",
                n_hidden, mass_min_da / 1000, mass_max_da / 1000)),
      tags$div(rows)
    )
  })

  # Recomputes on ANY iso_chk_*/iso_ptm_* change -- Shiny tracks dynamic
  # input[[...]] reads at runtime, so this reactive correctly re-fires even
  # though the input IDs themselves are generated dynamically above.
  derived <- reactive({
    cat_result <- catalog()
    pfs <- bare_pfs()
    # Returns an EMPTY result (not a hard req() stop) when nothing is loaded
    # yet -- e.g. right after the reset button clears catalog()/bare_pfs()
    # -- so outputs that depend on this (proteoform_table, digestion_result,
    # etc.) actually receive a fresh "nothing here" value and visibly clear
    # themselves, rather than req()'s default behavior of leaving whatever
    # they last rendered frozen on screen.
    if (is.null(cat_result) || is.null(pfs)) return(list(rows = list(), warnings = character(0)))
    rows <- list()
    warnings_all <- character(0)
    for (tid in cat_result$transcript_ids) {
      chk <- input[[paste0("iso_chk_", tid)]]
      if (is.null(chk) || !isTRUE(chk)) next
      pf <- pfs[[tid]]
      if (is.null(pf)) next
      bare_id <- paste0(tid, "#bare")
      rows[[bare_id]] <- list(id = bare_id, pf = pf, label = paste0(tid, " (unmodified)"), iso_key = tid)

      ptm_txt <- input[[paste0("iso_ptm_", tid)]]
      parsed <- parse_ptm_spec_text(pf$sequence, ptm_txt %||% "", tid)
      warnings_all <- c(warnings_all, parsed$warnings)
      for (gi in seq_along(parsed$groups)) {
        ptms <- parsed$groups[[gi]]
        pf_id <- paste0(tid, "#", gi)
        pf_obj <- proteoform(id = pf_id, sequence = pf$sequence, ptms = ptms, provenance = "module1_isoform_selection")
        label <- paste0(tid, " + ", paste(vapply(ptms, function(p) paste0(p$name, "@", p$site), character(1)), collapse = ", "))
        rows[[pf_id]] <- list(id = pf_id, pf = pf_obj, label = label, iso_key = tid)
      }
    }
    list(rows = rows, warnings = warnings_all)
  })

  output$ptm_warnings_ui <- renderUI({
    d <- derived()
    if (length(d$warnings) == 0) return(NULL)
    tags$div(lapply(d$warnings, function(w) tags$div(class = "pt-warn", w)))
  })

  output$proteoform_table <- renderTable({
    d <- derived()
    if (length(d$rows) == 0) return(data.frame(Proteoform = character(0), `Mass (Da)` = numeric(0), check.names = FALSE))
    data.frame(
      Proteoform = vapply(d$rows, function(r) r$label, character(1)),
      `Mass (Da)` = vapply(d$rows, function(r) round(proteoform_mass(r$pf)$mass, 2), numeric(1)),
      check.names = FALSE
    )
  })

  # ============================================================
  # Middle-down: simulate a limited-digestion protease treatment on whatever
  # full-length proteoforms Module 1/2/3 currently produced (derived()$rows),
  # sitting ABOVE those modules rather than belonging to any one of them --
  # same reasoning as the shared "MS resolution parameters" panel, since
  # "MS strategy" is selected even before that.
  # ============================================================
  digestion_result <- reactive({
    req(identical(input$ms_strategy, "middledown"))
    d <- derived()
    # Same reasoning as derived()'s own empty-safety above: an EMPTY result
    # (not req()) when there's nothing checked yet, so digestion-dependent
    # outputs clear themselves after a reset instead of freezing on stale
    # candidates from before the reset.
    if (length(d$rows) == 0) {
      return(list(candidates = NULL, peptides = list(), parent_label = character(0),
                  parent_iso_key = character(0), parent_counts = NULL))
    }
    pf_list <- lapply(d$rows, function(r) r$pf)
    labels <- vapply(d$rows, function(r) r$label, character(1))
    iso_keys <- vapply(d$rows, function(r) r$iso_key, character(1))
    mass_min_da <- (input$md_mass_min_kda %||% 3) * 1000
    mass_max_da <- (input$md_mass_max_kda %||% 10) * 1000
    digest_proteoform_set(
      pf_list, labels, iso_keys,
      enzyme = input$md_protease %||% "Lys-C",
      mass_min_da = mass_min_da, mass_max_da = mass_max_da,
      mode = input$ms_mode %||% "denatured",
      r_ref = input$ms_r_ref %||% 120000, mz_ref = input$ms_mz_ref %||% 200
    )
  })

  output$digestion_candidates_ui <- renderUI({
    dg <- digestion_result()
    if (is.null(dg$candidates) || nrow(dg$candidates) == 0) {
      return(tags$p(class = "pt-note", "No candidate peptides fall in the current mass window for this protease -- try widening the window or a different enzyme."))
    }
    cands <- dg$candidates
    # Default-check the single top-ranked (digest_proteoform() already sorts
    # by the feasibility-first, informativeness-second framework) candidate
    # per parent, as a starting point the user can freely override.
    default_ids <- unlist(lapply(split(cands$id, cands$parent_id), function(ids) ids[1]))

    # A real <table> (not the flex "pt-isorow" divs used for the shorter
    # isoform rows elsewhere) -- with 8 columns of real data this needs
    # proper table column-width negotiation and text wrapping WITHIN a
    # cell, not a flex row that just breaks onto a second visual line once
    # the modal is narrower than the sum of the columns' min-widths.
    header <- tags$tr(
      tags$th(style = "width:26px;"), tags$th("Parent"), tags$th("Range"),
      tags$th("Missed cl."), tags$th("Mass (Da)"), tags$th("MS1 FWHM (Da)"),
      tags$th("MS2 propensity"), tags$th("PTM sites")
    )
    rows <- lapply(seq_len(nrow(cands)), function(r) {
      cid <- cands$id[r]
      prev_checked <- isolate(input[[paste0("pep_chk_", sanitize_html_id(cid))]])
      default_checked <- if (!is.null(prev_checked)) prev_checked else (cid %in% default_ids)
      is_intact <- isTRUE(cands$is_intact[r])
      tags$tr(
        style = if (is_intact) "background:#fff6e5;" else NULL,
        tags$td(checkboxInput(paste0("pep_chk_", sanitize_html_id(cid)), NULL, value = default_checked)),
        tags$td(dg$parent_label[[cid]] %||% cands$parent_id[r]),
        tags$td(if (is_intact) sprintf("%d-%d (%d aa, intact -- no digest fragment in window)", cands$start[r], cands$end[r], cands$length[r])
                 else sprintf("%d-%d (%d aa)", cands$start[r], cands$end[r], cands$length[r])),
        tags$td(if (is_intact) "--" else cands$missed_cleavages[r]),
        tags$td(sprintf("%.1f", cands$mass[r])),
        tags$td(sprintf("%.2f", cands$ms1_fwhm_da[r])),
        tags$td(sprintf("%.2f", cands$ms2_avg_propensity[r])),
        tags$td(cands$ptm_sites_covered[r])
      )
    })
    tags$table(class = "pt-pep-table", tags$thead(header), tags$tbody(rows))
  })

  output$digestion_coverage_text <- renderText({
    dg <- digestion_result()
    req(dg$candidates)
    d <- derived()
    selected_ids <- Filter(function(cid) isTRUE(input[[paste0("pep_chk_", sanitize_html_id(cid))]]), dg$candidates$id)
    if (length(selected_ids) == 0) return("No candidates currently selected.")
    parent_lengths <- as.list(vapply(names(d$rows), function(pid) nchar(d$rows[[pid]]$pf$sequence), integer(1)))
    summary_df <- digestion_coverage_summary(dg$candidates, selected_ids, parent_lengths)
    if (nrow(summary_df) == 0) return("")
    paste(sprintf("%s: %d/%d residues covered (%.1f%%) by selected peptide(s).",
                  summary_df$parent_id, summary_df$covered_residues, summary_df$total_residues, summary_df$pct),
          collapse = "  ")
  })

  # One-line always-visible summary in the main page -- the full candidate
  # table (often dozens to 100+ rows per protein) lives in a modal instead,
  # opened on demand, so the main page doesn't force the user to scroll past
  # a huge table to reach the actual MS1/MS2 analysis panels below.
  output$digestion_summary_text <- renderText({
    dg <- digestion_result()
    counts <- dg$parent_counts
    n_intact <- if (!is.null(counts)) sum(counts$is_intact_only) else 0L
    intact_note <- if (n_intact > 0) sprintf(" (%d proteoform(s) had no digest fragment in this window for this protease, so the intact protein is offered instead -- see \"Select peptides...\".)", n_intact) else ""
    if (is.null(dg$candidates) || nrow(dg$candidates) == 0) {
      return(paste0("No candidate peptides in the current mass window -- try widening it or a different enzyme.", intact_note))
    }
    n_parents <- length(unique(dg$candidates$parent_id))
    selected_ids <- Filter(function(cid) isTRUE(input[[paste0("pep_chk_", sanitize_html_id(cid))]]), dg$candidates$id)
    sprintf("%d candidate peptide(s) across %d proteoform(s); %d selected for analysis.%s",
            nrow(dg$candidates), n_parents, length(selected_ids), intact_note)
  })

  # Per-parent candidate-count breakdown shown inside the picker modal --
  # a proteoform with a sparse-cleavage-site enzyme (e.g. OmpT, which only
  # cuts rare dibasic K/R-K/R sites) can legitimately have no real digest
  # fragment in the mass window even though it's checked and included;
  # digest_proteoform() falls back to offering the intact protein itself in
  # that case (see its doc comment) so the proteoform still has something
  # selectable, but that fallback is worth calling out explicitly here
  # rather than looking like an ordinary single-candidate result.
  output$digestion_per_parent_ui <- renderUI({
    dg <- digestion_result()
    if (is.null(dg$parent_counts) || nrow(dg$parent_counts) == 0) return(NULL)
    pc <- dg$parent_counts
    tags$div(style = "font-size:12.5px;color:#555;padding:4px 0;",
      tags$strong("Candidates found per proteoform: "),
      paste(mapply(function(label, n, is_intact) {
        if (is_intact) sprintf("%s (0 digest fragments in window -- intact protein offered instead)", label)
        else sprintf("%s (%d)", label, n)
      }, pc$label, pc$n_candidates, pc$is_intact_only), collapse = "; ")
    )
  })

  observeEvent(input$btn_open_peptide_picker, {
    showModal(modalDialog(
      title = "Select middle-down peptide candidates",
      size = "xl", easyClose = TRUE,
      p(class = "pt-note", "Every in-silico digest fragment (any number of missed cleavages) of the checked proteoforms above that falls in the mass window. Sorted by likely feasibility first (fewer missed cleavages, tighter MS1 peak), then by how many PTM sites it covers. Pick which candidate(s) to treat as \"proteins\" for the MS1/MS2 analysis below -- the top-ranked candidate per parent is pre-checked as a starting point."),
      uiOutput("digestion_per_parent_ui"),
      div(style = "max-height:55vh;overflow:auto;", uiOutput("digestion_candidates_ui")),
      hr(),
      textOutput("digestion_coverage_text"),
      footer = modalButton("Done")
    ))
  })

  observe({
    if (identical(input$ms_strategy, "middledown")) {
      dg <- digestion_result()
      if (is.null(dg$candidates) || nrow(dg$candidates) == 0) {
        updateSelectInput(session, "pf_target_select", choices = character(0))
        return()
      }
      selected_ids <- Filter(function(cid) isTRUE(input[[paste0("pep_chk_", sanitize_html_id(cid))]]), dg$candidates$id)
      if (length(selected_ids) == 0) {
        updateSelectInput(session, "pf_target_select", choices = character(0))
        return()
      }
      labels <- vapply(selected_ids, function(cid) {
        row <- dg$candidates[dg$candidates$id == cid, ]
        sprintf("%s [%d-%d, %d aa]", dg$parent_label[[cid]] %||% cid, row$start[1], row$end[1], row$length[1])
      }, character(1))
      updateSelectInput(session, "pf_target_select", choices = setNames(selected_ids, labels))
    } else {
      d <- derived()
      choices <- if (length(d$rows) == 0) character(0) else {
        setNames(names(d$rows), vapply(d$rows, function(r) r$label, character(1)))
      }
      updateSelectInput(session, "pf_target_select", choices = choices)
    }
  })

  # "Run analysis" only makes sense once there's actually something in the
  # proteoform table AND a confounder-search target has been picked -- grey
  # it out otherwise so a click can't run on nothing or with no target.
  # Top-down: at least one checked isoform. Middle-down: at least one
  # checked peptide candidate (a checked proteoform with zero digest
  # candidates in the window doesn't count as "content" on its own, but see
  # digest_proteoform()'s intact-protein fallback -- that fallback always
  # gives it something checkable, so this still resolves TRUE for that case
  # once picked). Also requires both "Set" buttons above (MS strategy, MS
  # resolution parameters) to be locked in, same as Load isoforms/Run
  # sequence analysis, since every one of these entry points depends on
  # those global settings.
  observe({
    has_content <- if (identical(input$ms_strategy, "middledown")) {
      dg <- digestion_result()
      if (is.null(dg$candidates) || nrow(dg$candidates) == 0) {
        FALSE
      } else {
        length(Filter(function(cid) isTRUE(input[[paste0("pep_chk_", sanitize_html_id(cid))]]), dg$candidates$id)) > 0
      }
    } else {
      length(derived()$rows) > 0
    }
    has_target <- !is.null(input$pf_target_select) && nzchar(input$pf_target_select)
    enabled <- has_content && has_target && ms_strategy_set() && ms_resolution_set()
    session$sendCustomMessage("pt_set_button_enabled", list(id = "btn_run_analysis", enabled = enabled))
  })

  # analysis_locked/analysis_ready MUST be set from an observeEvent of their
  # own, not from inside the `analysis` eventReactive body below: an
  # eventReactive only actually evaluates when something reads its value,
  # and viz_script/confounder_status's own analysis_ready() guards (below)
  # mean they may never call analysis() at all while analysis_ready() is
  # still FALSE -- setting it TRUE only inside analysis()'s body would then
  # be a deadlock (analysis() never runs because nothing reads it while
  # analysis_ready() is FALSE, and analysis_ready() never becomes TRUE
  # because analysis() never runs). This observeEvent fires unconditionally
  # on every click instead, independent of that.
  observeEvent(input$btn_run_analysis, {
    analysis_locked(TRUE)
    analysis_ready(TRUE)
  })

  # ---- Run-gated heavy analysis (MS1 + ladder + confounder search) ----
  analysis <- eventReactive(input$btn_run_analysis, {
    r_ref <- isolate(input$ms_r_ref)
    mz_ref <- isolate(input$ms_mz_ref)
    safety_margin <- isolate(input$ms_safety_margin)
    mode <- isolate(input$ms_mode)
    strategy <- isolate(input$ms_strategy)
    protease <- isolate(input$md_protease)
    cat_result <- isolate(catalog())
    residue_offset <- list()
    dg <- NULL

    # Real confounders for a middle-down target are OTHER proteins' digest
    # peptides (cut with the SAME protease), not other intact proteins -- an
    # intact 4 kDa protein and a 4 kDa peptide cut from the middle of some
    # unrelated 60 kDa protein are both real collision risks, but only a
    # digested pool captures the second. Falls back to the intact-protein
    # index for top-down. First use of a given enzyme is a real ~20-40s
    # proteome-scale computation (cached to disk after that), hence the
    # progress message -- see get_digested_reference_pool().
    confounder_mass_index <- if (identical(strategy, "middledown")) {
      withProgress(message = paste("Preparing", protease, "-digested confounder pool..."), value = 0.3, {
        tryCatch(
          get_digested_reference_pool(protease, on_build_start = function() {
            incProgress(0, detail = "first use of this enzyme: one-time proteome-wide digest (~30s), then cached to disk")
          }),
          error = function(e) NULL
        )
      })
    } else {
      # Top-down mass filter (see "MS strategy" panel): the confounder
      # search pool is restricted to the SAME intact-protein mass range the
      # isoform catalog above is already limited to, so a candidate outside
      # the instrument's usable mass range can't turn up as a "confounder"
      # the user could never actually observe anyway.
      topdown_mass_min_da <- (isolate(input$topdown_mass_min_kda) %||% 10) * 1000
      topdown_mass_max_da <- (isolate(input$topdown_mass_max_kda) %||% 220) * 1000
      reference_mass_index[reference_mass_index$mass >= topdown_mass_min_da &
                              reference_mass_index$mass <= topdown_mass_max_da, ]
    }

    if (identical(strategy, "middledown")) {
      dg <- isolate(digestion_result())
      req(dg$candidates, nrow(dg$candidates) > 0)
      selected_ids <- Filter(function(cid) isTRUE(isolate(input[[paste0("pep_chk_", sanitize_html_id(cid))]])), dg$candidates$id)
      req(length(selected_ids) > 0)
      pf_list <- dg$peptides[selected_ids]
      iso_key_of <- setNames(dg$parent_iso_key[selected_ids], selected_ids)
      for (cid in selected_ids) {
        row <- dg$candidates[dg$candidates$id == cid, ]
        residue_offset[[cid]] <- row$start[1] - 1L
      }
    } else {
      d <- isolate(derived())
      checked_ids <- names(d$rows)  # everything currently in the proteoform table is "included"
      req(length(checked_ids) > 0)
      pf_list <- lapply(d$rows[checked_ids], function(r) r$pf)
      names(pf_list) <- checked_ids
      iso_key_of <- vapply(d$rows[checked_ids], function(r) r$iso_key, character(1))
    }

    masses <- vapply(pf_list, function(pf) proteoform_mass(pf)$mass, numeric(1))
    tiers <- compute_ladder_tiers(pf_list, r_ref = r_ref, mz_ref = mz_ref, safety_margin = safety_margin, mode = mode)

    payload1 <- build_section1_payload(pf_list, iso_key_of, masses, tiers, cat_result$exon_table, residue_offset = residue_offset)
    for (i in seq_along(payload1$proteoforms)) {
      pf <- pf_list[[payload1$proteoforms[[i]]$id]]
      env <- predict_charge_envelope(pf$sequence, masses[[payload1$proteoforms[[i]]$id]], mode = mode)
      payload1$proteoforms[[i]]$env <- Map(function(z, mz, ri) list(z = z, mz = round(mz, 2), rel = round(ri, 4)),
                                            env$z, env$mz, env$relative_intensity)
    }

    target_id <- isolate(input$pf_target_select)
    payload2 <- NULL
    if (!is.null(target_id) && nzchar(target_id) && target_id %in% names(pf_list) && !is.null(confounder_mass_index)) {
      target_pf <- pf_list[[target_id]]
      target_mass <- masses[[target_id]]
      confounder_result <- tryCatch(
        search_confounding_proteins_real(target_pf, confounder_mass_index, mode = mode,
                                          r_ref = r_ref, mz_ref = mz_ref, safety_margin = safety_margin),
        error = function(e) NULL
      )
      target_iso_key <- iso_key_of[[target_id]]
      target_exon_table <- if (!is.null(cat_result)) cat_result$exon_table[cat_result$exon_table$transcript_id == target_iso_key, ] else NULL
      if (identical(strategy, "middledown") && !is.null(target_exon_table) && nrow(target_exon_table) > 0) {
        row <- dg$candidates[dg$candidates$id == target_id, ]
        target_exon_table <- shift_exon_table_for_peptide(target_exon_table, row$start[1], row$end[1])
      }
      confounder_tiers <- compute_confounder_tiers(target_pf, list(), r_ref = r_ref, mz_ref = mz_ref, safety_margin = safety_margin, mode = mode)
      confounder_envs <- list()
      if (!is.null(confounder_result) && nrow(confounder_result$candidates) > 0) {
        conf_pfs <- build_confounder_proteoforms(confounder_result$candidates$id, confounder_mass_index)
        confounder_tiers <- compute_confounder_tiers(target_pf, conf_pfs, r_ref = r_ref, mz_ref = mz_ref, safety_margin = safety_margin, mode = mode)
        confounder_envs <- lapply(conf_pfs, function(p) predict_charge_envelope(p$sequence, proteoform_mass(p)$mass, mode = mode))
      }
      payload2 <- build_section2_payload(target_pf, target_mass, confounder_tiers, target_exon_table, confounder_result, confounder_envs)
      target_env <- predict_charge_envelope(target_pf$sequence, target_mass, mode = mode)
      payload2$target$env <- Map(function(z, mz, ri) list(z = z, mz = round(mz, 2), rel = round(ri, 4)),
                                  target_env$z, target_env$mz, target_env$relative_intensity)
    }

    list(payload1 = payload1, payload2 = payload2, target_id = target_id)
  })

  output$stale_notice <- renderText({
    # Empty-safe (not req()) so the reset button's catalog(NULL) actually
    # clears this text instead of req()'s default "silently keep whatever
    # was last shown" -- same reasoning as derived()'s empty-safety above.
    if (is.null(catalog())) return("")
    "Click \"Run analysis\" to (re)compute MS1/MS2/confounder results for the currently included proteoforms."
  })

  output$viz_script <- renderUI({
    # Guard on catalog() (a plain reactiveVal) BEFORE touching analysis()
    # (an eventReactive, frozen at its last click) -- reset sets catalog()
    # to NULL, and this stops a stale script tag (from before the reset)
    # from being reinserted, without needing to overwrite this reactive
    # binding itself (which would break it for every analysis after the
    # first reset -- eventReactive() has no way to be told "return empty").
    # analysis_ready() additionally guards against replaying analysis()'s
    # stale frozen value for data loaded AFTER a reset but before "Run
    # analysis" is clicked again (see its own doc comment above).
    if (is.null(catalog()) || !analysis_ready()) return(NULL)
    a <- analysis()
    j1 <- jsonlite::toJSON(a$payload1, auto_unbox = TRUE, digits = 4, null = "null")
    j2 <- if (!is.null(a$payload2)) jsonlite::toJSON(a$payload2, auto_unbox = TRUE, digits = 4, null = "null") else "null"
    tags$script(HTML(sprintf(
      "PT.renderSection1(%s); var __pt_p2 = %s; if (__pt_p2) { PT.renderSection2(__pt_p2); } else { document.getElementById('s2-ladder').innerHTML=''; document.getElementById('s2-ms1').innerHTML=''; document.getElementById('s2-legend').innerHTML='<p style=\"color:#888;font-size:12px;\">No confounder-search target selected, or no reference mass index loaded.</p>'; }",
      j1, j2
    )))
  })

  output$confounder_status <- renderText({
    if (is.null(catalog()) || !analysis_ready()) return("")
    a <- analysis()
    if (is.null(a$target_id) || !nzchar(a$target_id)) return("Pick a confounder-search target below the proteoform table, then Run analysis.")
    if (is.null(reference_mass_index)) return("No reference-proteome mass index loaded -- run scripts/build_reference_proteome.R.")
    if (is.null(a$payload2)) return("Confounder search failed for this target.")
    sprintf("Target: %s. %d real confounding protein(s) found in window.", a$target_id, length(a$payload2$confounders))
  })

  # ============================================================
  # Option 2: FASTA upload. Two genuinely independent steps:
  #   1. run_fasta_alignment(): minimap2 spliced alignment -> which known
  #      gene/isoforms this sequence belongs to. Never looks at translation.
  #   2. run_fasta_translation(): TransDecoder's top-scoring ORF candidates.
  #      Never looks at the genome alignment.
  # The user then picks a known isoform to compare exon structure against,
  # picks which ORF candidate is the real translation, adds PTMs, and only
  # THEN (via "Add to comparison") does any of this feed into the SAME
  # catalog()/bare_pfs() reactives Option 1 uses -- so the rest of the app
  # (proteoform table, Run analysis, MS1/MS2, confounder search) works
  # unchanged for a novel sequence too.
  # ============================================================
  fasta_alignment <- reactiveVal(NULL)
  fasta_orfs <- reactiveVal(NULL)

  # Same "paste box or uploaded file, whichever is present" precedence the
  # actual run handler below uses, but read reactively here purely to gate
  # "Run sequence analysis" -- greyed out until there's real FASTA content
  # (not just an empty/garbage paste or an unselected file input).
  current_fasta_text <- reactive({
    if (!is.null(input$fasta_file)) {
      tryCatch(paste(readLines(input$fasta_file$datapath), collapse = "\n"), error = function(e) "")
    } else {
      input$fasta_text %||% ""
    }
  })
  observe({
    session$sendCustomMessage("pt_set_button_enabled", list(
      id = "btn_run_fasta",
      enabled = looks_like_valid_fasta(current_fasta_text()) && ms_strategy_set() && ms_resolution_set()
    ))
  })

  observeEvent(input$btn_run_fasta, {
    analysis_locked(TRUE)
    mmi_path <- "reference/genome/GRCh38.mmi"
    if (!file.exists(mmi_path) || file.size(mmi_path) == 0) {
      output$fasta_status <- renderText(
        "Reference genome index is still building in the background (reference/genome/GRCh38.mmi). This pathway will activate automatically once it's ready."
      )
      return()
    }

    fasta_text <- if (!is.null(input$fasta_file)) {
      paste(readLines(input$fasta_file$datapath), collapse = "\n")
    } else {
      input$fasta_text
    }
    if (is.null(fasta_text) || !nzchar(trimws(fasta_text %||% ""))) {
      output$fasta_status <- renderText("Paste a FASTA sequence or upload a file first.")
      return()
    }

    withProgress(message = "Aligning (minimap2) and translating (TransDecoder)...", value = 0.2, {
      aln <- tryCatch(run_fasta_alignment(fasta_text, reference_exon_index, mmi_path = mmi_path, transcript_id = "NOVEL_1"),
                       error = function(e) list(align_error = conditionMessage(e)))
      incProgress(0.5)
      orfs <- tryCatch(run_fasta_translation(fasta_text, transcript_id = "NOVEL_1"),
                        error = function(e) list(orf_error = conditionMessage(e)))
      fasta_alignment(aln)
      fasta_orfs(orfs)
    })

    a <- fasta_alignment()
    if (!is.null(a$align_error)) {
      output$fasta_status <- renderText(paste("minimap2 alignment failed:", a$align_error))
      return()
    }
    status <- if (!is.null(a$matched_gene)) {
      sprintf("minimap2: matched to %s (%d overlapping exons, %d known transcripts) -- pick one below to compare against.",
              a$matched_gene$gene_name, a$matched_gene$n_overlap_exons, length(a$matched_gene$matched_transcript_ids))
    } else {
      paste("minimap2:", paste(a$warnings, collapse = " "))
    }
    o <- fasta_orfs()
    if (!is.null(o$orf_error)) {
      status <- paste(status, "TransDecoder failed:", o$orf_error)
    } else if (nrow(o) == 0) {
      status <- paste(status, "TransDecoder found no candidate ORFs.")
    } else {
      status <- paste(status, sprintf("TransDecoder: %d candidate ORF(s) found, top-scoring shown first.", nrow(o)))
    }
    output$fasta_status <- renderText(status)
  })

  output$fasta_results_ui <- renderUI({
    a <- fasta_alignment()
    # Empty-safe (not req()) so the reset button's fasta_alignment(NULL)
    # actually clears this panel instead of freezing on stale ORF/exon
    # content from before the reset.
    if (is.null(a) || !is.null(a$align_error)) return(NULL)

    isoform_choices <- character(0)
    if (!is.null(a$matched_gene)) {
      known_cat <- build_gene_isoform_catalog(a$matched_gene$gene_name, reference_exon_index)
      isoform_choices <- setNames(known_cat$transcript_ids,
                                    sprintf("%s (%d aa)", known_cat$transcript_ids, known_cat$protein_lengths))
    }

    o <- fasta_orfs()
    orf_choices <- character(0)
    if (!is.null(o) && is.null(o$orf_error) && nrow(o) > 0) {
      orf_choices <- setNames(as.character(seq_len(nrow(o))),
                                sprintf("%s -- %s, %d aa, score %.1f", o$orf_id, o$quality, o$aa_len, o$score))
    }

    tagList(
      hr(),
      h5("1. TransDecoder candidate translations"),
      p(class = "pt-note", "Pick which candidate is the real translation -- this decides the coding (CDS) vs. non-coding (UTR) split shown for your sequence in the exon alignment below."),
      tableOutput("fasta_orf_table"),
      if (length(orf_choices) > 0) tagList(
        selectInput("fasta_orf_select", "Use this translation", choices = orf_choices, width = "420px"),
        textInput("fasta_ptm_spec", "PTMs on this translation (optional)",
                  placeholder = "e.g. 133_Thr_Phospho; 210_Pro_Oxidation,215_Ser_Sulfo", width = "420px")
      ) else tags$p(class = "pt-note", "No ORF candidates to select."),
      hr(),
      h5("2. Known isoform(s) to compare against"),
      if (length(isoform_choices) > 0) tagList(
        p(class = "pt-note", "Select one or more known transcripts of the matched gene (ctrl/cmd-click, or use the dropdown, to select multiple)."),
        selectInput("fasta_isoform_select", "Compare against", choices = isoform_choices,
                     selected = isoform_choices[1], multiple = TRUE, width = "420px")
      ) else tags$p(class = "pt-note", "No known gene matched -- no known isoform to compare exon structure against."),
      hr(),
      h5("3. Exon structure alignment"),
      p(class = "pt-note", "Reflects your choices above: coding (CDS) vs. non-coding (UTR) for your sequence follows the translation picked in step 1; the rows compared follow step 2. Scroll/pinch or use the +/- buttons to zoom, drag to pan."),
      div(id = "fasta-exon-zoom"),
      tags$svg(id = "fasta-exon-align", class = "pt-viz", viewBox = "0 0 640 100"),
      div(id = "fasta-exon-align-legend"),
      hr(),
      if (length(orf_choices) > 0) actionButton("btn_fasta_add_to_comparison", "Add to comparison", class = "btn-success")
    )
  })

  output$fasta_orf_table <- renderTable({
    o <- fasta_orfs()
    req(o, is.null(o$orf_error))
    if (nrow(o) == 0) return(data.frame(message = "No candidate ORFs found."))
    data.frame(
      Candidate = o$orf_id, `Start codon` = ifelse(o$has_start, "Yes", "No"),
      `Stop codon` = ifelse(o$has_stop, "Yes", "No"), `AA count` = o$aa_len,
      Quality = o$quality, Score = round(o$score, 1),
      check.names = FALSE
    )
  })

  observeEvent(list(input$fasta_isoform_select, input$fasta_orf_select), {
    a <- fasta_alignment()
    req(a, is.null(a$align_error))
    preview <- if (length(input$fasta_isoform_select) == 0) {
      NULL
    } else {
      # The novel row's own coding/non-coding split depends on whichever
      # TransDecoder candidate is currently selected -- gene-matching and
      # translation stay independent (see run_fasta_alignment()'s doc
      # comment), but marking coding sequence is inherently a translation
      # question, so it's the one place this view reads the ORF selection.
      novel_cds <- NULL
      o <- fasta_orfs()
      if (!is.null(o) && is.null(o$orf_error) && nrow(o) > 0 && !is.null(input$fasta_orf_select)) {
        orf_row <- o[as.integer(input$fasta_orf_select), ]
        gene_name <- if (!is.null(a$matched_gene)) a$matched_gene$gene_name else NA_character_
        novel_cds <- tryCatch(build_selected_orf_exon_table(orf_row, a, "NOVEL_1", gene_name), error = function(e) NULL)
      }
      tryCatch(
        build_exon_alignment_preview(a, input$fasta_isoform_select, reference_exon_index, novel_cds_range = novel_cds),
        error = function(e) NULL
      )
    }
    session$sendCustomMessage("pt_render_exon_alignment", list(has_data = !is.null(preview), payload = preview))
  }, ignoreNULL = FALSE)

  observeEvent(input$btn_fasta_add_to_comparison, {
    a <- fasta_alignment()
    o <- fasta_orfs()
    req(a, o, is.null(a$align_error), is.null(o$orf_error), input$fasta_orf_select)
    orf_row <- o[as.integer(input$fasta_orf_select), ]

    gene_name <- if (!is.null(a$matched_gene)) a$matched_gene$gene_name else NA_character_
    novel_exon_table <- build_selected_orf_exon_table(orf_row, a, "NOVEL_1", gene_name)
    empty_exon_cols <- c("transcript_id", "gene_id", "gene_name", "is_canonical", "exon_number",
                          "seqname", "strand", "start", "end", "residue_start", "residue_end", "protein_length")
    if (is.null(novel_exon_table)) {
      novel_exon_table <- setNames(data.frame(matrix(nrow = 0, ncol = length(empty_exon_cols))), empty_exon_cols)
    }

    novel_pf <- proteoform(id = "NOVEL_1", sequence = orf_row$protein_sequence, provenance = "module2_longread_orf")

    # NOVEL_1's own PTM spec textbox (rendered by isoform_catalog_ui, same
    # as any Option 1 transcript) will parse whatever's prefilled here into
    # PTM-variant proteoforms via the existing derived() reactive -- no
    # need to build them by hand here.
    ptm_prefill(list(NOVEL_1 = input$fasta_ptm_spec %||% ""))

    extra_checked <- "NOVEL_1"
    if (length(input$fasta_isoform_select) > 0) {
      selected_tids <- input$fasta_isoform_select
      # minimap2/exon alignment above stays transcript-level (that's
      # inherently genomic/splicing structure), but converting the selected
      # canonical transcripts into "proteins" for the MS1/MS2 analysis below
      # collapses any that share an identical sequence to one representative
      # -- same reasoning/helper as Module 1's isoform catalog.
      selected_pf_raw <- build_proteoforms_for_transcripts(selected_tids)
      deduped <- dedupe_proteoforms_by_sequence(selected_pf_raw)
      selected_pf <- deduped$pf_list
      isoform_synonyms(deduped$synonyms)
      selected_exon_table <- reference_exon_index[reference_exon_index$transcript_id %in% selected_tids, ]

      merged_exon_table <- rbind(novel_exon_table, selected_exon_table)
      merged_transcript_ids <- c("NOVEL_1", names(selected_pf))
      merged_protein_lengths <- c(NOVEL_1 = nchar(orf_row$protein_sequence),
                                    setNames(vapply(selected_pf, function(p) nchar(p$sequence), integer(1)), names(selected_pf)))

      catalog(list(gene = gene_name, exon_table = merged_exon_table,
                    transcript_ids = merged_transcript_ids, protein_lengths = merged_protein_lengths))
      bare_pfs(c(setNames(list(novel_pf), "NOVEL_1"), selected_pf))
      extra_checked <- c(extra_checked, names(selected_pf))
    } else {
      isoform_synonyms(list())
      catalog(list(gene = gene_name, exon_table = novel_exon_table,
                    transcript_ids = "NOVEL_1", protein_lengths = c(NOVEL_1 = nchar(orf_row$protein_sequence))))
      bare_pfs(setNames(list(novel_pf), "NOVEL_1"))
    }
    precheck_extra(extra_checked)
    catalog_source("fasta")

    output$gene_status <- renderText(sprintf(
      "From FASTA input: %s translation (%d aa) added%s. Pick it (and any known isoforms) below, then Run analysis.",
      orf_row$orf_id, nchar(orf_row$protein_sequence),
      if (!is.na(gene_name)) paste0(" alongside ", gene_name, "'s known isoforms") else ""
    ))
    # Reveals the shared proteoform/MS1/MS2 panel WITHOUT switching
    # input_mode itself -- flipping input_mode to "gene" here used to make
    # the visible radio jump to "1. Gene -> isoform -> proteoform", which
    # reads as though the app had left module 2, even though the data on
    # screen is still the FASTA-derived one. force_gene_view (a hidden
    # checkbox, see ui.R) drives the SAME panel's conditionalPanel instead,
    # so the radio keeps showing "2. FASTA sequence" as selected.
    updateCheckboxInput(session, "force_gene_view", value = TRUE)
  })

  # ============================================================
  # Option 3: rMATS upload. rMATS only reports the differential exon(s) and
  # their immediate flanking exons, not the rest of the transcript -- so
  # instead of trying to invent a full transcript from a handful of exons'
  # worth of information, match_rmats_se_transcripts()/
  # match_rmats_mxe_transcripts() (R/rmats_adapter.R) look up which
  # already-annotated transcripts of the gene (same reference_exon_index
  # Option 1 uses) structurally match each "arm" of the event (e.g.
  # exon-inclusion vs exon-skipping for SE; 1st-exon vs 2nd-exon for MXE).
  # The user picks which of those real transcripts to send into the SAME
  # catalog()/bare_pfs() pathway Option 1/2 use -- so the rest of the app
  # (proteoform table, digestion, Run analysis, MS1/MS2, confounder search)
  # works unchanged here too. Only SE and MXE are supported so far.
  # ============================================================
  rmats_event <- reactiveVal(NULL)
  rmats_matches <- reactiveVal(NULL)

  # Same "there's a real file selected" gate as Option 2's content check --
  # rMATS files are always uploaded (no paste option), and actual
  # parseability is cheap enough to just check at click time (see the run
  # handler below) rather than duplicating the parser here for gating.
  observe({
    session$sendCustomMessage("pt_set_button_enabled", list(
      id = "btn_run_rmats",
      enabled = !is.null(input$rmats_file) && input$rmats_event_type %in% c("SE", "MXE") &&
        ms_strategy_set() && ms_resolution_set()
    ))
  })

  observeEvent(input$btn_run_rmats, {
    event_type <- input$rmats_event_type
    req(input$rmats_file, event_type %in% c("SE", "MXE"))
    parser <- if (identical(event_type, "MXE")) parse_rmats_mxe else parse_rmats_se
    matcher <- if (identical(event_type, "MXE")) match_rmats_mxe_transcripts else match_rmats_se_transcripts
    events <- tryCatch(parser(input$rmats_file$datapath), error = function(e) NULL)
    if (is.null(events) || nrow(events) == 0) {
      output$rmats_status <- renderText(sprintf(
        "Could not parse this file as rMATS %s results -- check it's the %s.MATS.JC/JCEC output with its header row intact.", event_type, event_type
      ))
      rmats_event(NULL)
      rmats_matches(NULL)
      return()
    }
    event <- events[1, ]
    rmats_event(event)
    matches <- matcher(event, reference_exon_index)
    rmats_matches(matches)

    arm_counts_text <- paste(vapply(matches$arms, function(a) sprintf("%d match the %s", nrow(a$candidates), a$label), character(1)), collapse = "; ")
    multi_note <- if (nrow(events) > 1) sprintf(" (file has %d events; only the first is used per run for now)", nrow(events)) else ""
    output$rmats_status <- renderText(sprintf(
      "%s %s event, %s (%s strand)%s: %s. Pick which to include below, then Add to comparison.",
      event$gene_symbol, event_type, rmats_event_region_text(event, event_type), event$strand, multi_note, arm_counts_text
    ))

    all_ids <- unique(unlist(lapply(matches$arms, function(a) a$candidates$transcript_id)))
    if (length(all_ids) > 0) {
      preview <- build_rmats_candidate_alignment(matches, reference_exon_index)
      session$sendCustomMessage("pt_render_exon_alignment", list(
        has_data = !is.null(preview), payload = preview, instance = "rmats",
        svg_id = "rmats-exon-align", legend_id = "rmats-exon-align-legend",
        zoom_id = "rmats-exon-zoom", id_prefix = "rmats-exon"
      ))
    } else {
      session$sendCustomMessage("pt_render_exon_alignment", list(has_data = FALSE, instance = "rmats",
        svg_id = "rmats-exon-align", legend_id = "rmats-exon-align-legend", zoom_id = "rmats-exon-zoom"))
    }
  })

  output$rmats_results_ui <- renderUI({
    matches <- rmats_matches()
    event <- rmats_event()
    if (is.null(event) || is.null(matches)) return(NULL)

    # Pre-check the single best-matching pair (highest pairing_score) per
    # arm as a sensible default rather than an arbitrary one -- see
    # match_rmats_arm_transcripts()'s doc comment for why pairing_score
    # (exon-set similarity against candidates in OTHER arms) is the right
    # criterion: it favors the pair that most likely represents "the same
    # underlying transcript, +/- this event" over transcripts that only
    # coincidentally match locally.
    default_checked <- unlist(lapply(matches$arms, function(a) {
      if (nrow(a$candidates) > 0) a$candidates$transcript_id[1] else NULL
    }))

    row_for <- function(tid, is_canonical, pairing_score, arm_label) {
      prev_checked <- isolate(input[[paste0("rmats_chk_", tid)]])
      checked <- if (!is.null(prev_checked)) prev_checked else (tid %in% default_checked)
      tags$div(class = "pt-isorow",
        checkboxInput(paste0("rmats_chk_", tid), NULL, value = checked),
        tags$span(class = "pt-id", tid),
        tags$span(class = "pt-meta", if (is_canonical) "canonical" else ""),
        tags$span(class = "pt-meta", sprintf("pairing score %.2f", pairing_score)),
        tags$span(class = "pt-note", sprintf("(%s)", arm_label))
      )
    }
    arm_blocks <- lapply(matches$arms, function(a) {
      rows <- if (nrow(a$candidates) > 0) {
        lapply(seq_len(nrow(a$candidates)), function(i) {
          row_for(a$candidates$transcript_id[i], a$candidates$is_canonical[i], a$candidates$pairing_score[i], a$label)
        })
      } else list()
      tagList(
        h6(a$label),
        if (length(rows) > 0) tags$div(rows) else tags$p(class = "pt-note", sprintf("No annotated transcript matches the %s.", tolower(a$label)))
      )
    })
    total_candidates <- sum(vapply(matches$arms, function(a) nrow(a$candidates), integer(1)))

    tagList(
      hr(),
      h5("1. Matching transcripts"),
      p(class = "pt-note", "Pairing score = how much of each transcript's exon structure (outside this event) agrees with its best match in another arm -- close to 1 means \"the same underlying transcript, using this arm's exon(s) vs. another's\"; the best-scoring pair is pre-checked below."),
      arm_blocks,
      hr(),
      h5("2. Exon structure alignment"),
      p(class = "pt-note", "Every matched transcript shown for context, regardless of which are checked above; the event's own differential exon(s) are boxed on top of the usual common/partial/unique coloring. Scroll/pinch or use the +/- buttons to zoom, drag to pan. (Rendered below, outside this dynamic panel.)"),
      hr(),
      if (total_candidates > 0) actionButton("btn_rmats_add_to_comparison", "Add to comparison", class = "btn-success")
    )
  })

  observeEvent(input$btn_rmats_add_to_comparison, {
    matches <- rmats_matches()
    req(matches)
    all_ids <- unique(unlist(lapply(matches$arms, function(a) a$candidates$transcript_id)))
    checked_ids <- Filter(function(tid) isTRUE(input[[paste0("rmats_chk_", tid)]]), all_ids)
    req(length(checked_ids) > 0)

    selected_pf_raw <- build_proteoforms_for_transcripts(checked_ids)
    req(length(selected_pf_raw) > 0)
    # Multiple checked transcripts commonly translate to the IDENTICAL
    # protein (same reasoning as Option 1's own isoform catalog -- see
    # dedupe_proteoforms_by_sequence()'s doc comment) -- collapsed here too,
    # so "N checked" and "N rows in the proteoform table" can legitimately
    # differ; the status message below says so explicitly (mirroring Option
    # 1's own "collapsed to N unique protein sequence(s)" wording) so this
    # doesn't read as checked selections silently going missing.
    deduped <- dedupe_proteoforms_by_sequence(selected_pf_raw)
    isoform_synonyms(deduped$synonyms)

    event <- rmats_event()
    event_type <- input$rmats_event_type
    exon_table <- reference_exon_index[reference_exon_index$transcript_id %in% names(deduped$pf_list), ]
    protein_lengths <- vapply(deduped$pf_list, function(p) nchar(p$sequence), integer(1))

    catalog(list(gene = event$gene_symbol, exon_table = exon_table,
                  transcript_ids = names(deduped$pf_list), protein_lengths = protein_lengths))
    bare_pfs(deduped$pf_list)
    catalog_source("rmats")
    precheck_extra(names(deduped$pf_list))

    n_checked <- length(checked_ids)
    n_unique <- length(deduped$pf_list)
    output$gene_status <- renderText(sprintf(
      "From rMATS %s event (%s, %s): %d transcript(s) selected, %s. Pick which to include below, then Run analysis.",
      event_type, event$gene_symbol, rmats_event_region_text(event, event_type), n_checked,
      if (n_unique < n_checked) sprintf("collapsed to %d unique protein sequence(s)", n_unique) else sprintf("%d unique protein sequence(s)", n_unique)
    ))
    # Same reasoning as Option 2's "Add to comparison" -- reveal the shared
    # proteoform/MS1/MS2 panel without switching input_mode itself, so the
    # radio keeps showing "3. rMATS..." as selected.
    updateCheckboxInput(session, "force_gene_view", value = TRUE)
  })

  # ============================================================
  # Reset: clear every module's input/loaded state back to a blank page, so
  # starting a new analysis doesn't require manually clearing text boxes and
  # unchecking isoforms one at a time. Deliberately does NOT touch "MS
  # strategy"/"MS resolution parameters" -- those are global settings, not
  # per-analysis state, same distinction the app already draws elsewhere.
  # ============================================================
  observeEvent(input$btn_reset_all, {
    # Module 1 (gene -> isoform -> proteoform). Nulling catalog()/bare_pfs()
    # cascades proteoform_table/isoform_catalog_ui/digestion_*/stale_notice/
    # viz_script/confounder_status back to their empty states on its own now
    # (see derived()/digestion_result() and those three outputs' own
    # empty-safety checks above) -- gene_status is the only one of this
    # group assigned directly (inside observeEvent blocks, not a reactive
    # expression tied to catalog()), so it's the only one that needs
    # clearing here too. Do NOT reassign stale_notice/viz_script/
    # confounder_status/fasta_results_ui here: those ARE reactive
    # expressions (renderText({...})/renderUI({...}) bound to catalog()/
    # analysis()/fasta_alignment()), and overwriting the binding itself
    # would permanently break them for the rest of the session -- the next
    # gene loaded or "Run analysis" click would have nothing to re-render.
    catalog(NULL)
    bare_pfs(NULL)
    precheck_extra(character(0))
    ptm_prefill(list())
    isoform_synonyms(list())
    catalog_source(NULL)
    analysis_locked(FALSE)
    analysis_ready(FALSE)
    updateTextInput(session, "gene_symbol", value = "")
    output$gene_status <- renderText("")
    updateSelectInput(session, "pf_target_select", choices = character(0))
    updateCheckboxInput(session, "force_gene_view", value = FALSE)

    # Module 2 (FASTA)
    fasta_alignment(NULL)
    fasta_orfs(NULL)
    updateTextAreaInput(session, "fasta_text", value = "")
    output$fasta_status <- renderText("")

    # Module 3 (rMATS)
    rmats_event(NULL)
    rmats_matches(NULL)
    output$rmats_status <- renderText("")

    # Raw SVG/legend/filter/zoom DOM content a PREVIOUS script-injected
    # render left behind -- clearing the render functions above stops
    # FUTURE updates but does not undo markup an already-executed <script>
    # tag wrote directly into these elements' innerHTML. Payload must be
    # non-empty -- jsonlite::toJSON(list()) serializes to "[]", and a
    # message sent with that empty-list payload silently never reached the
    # client at all (confirmed directly: server-side logging showed
    # sendCustomMessage() was reached and returned normally, but the
    # browser-side handler never fired) -- some part of Shiny's own custom
    # message dispatch treats a completely empty payload as nothing to send.
    session$sendCustomMessage("pt_reset_analysis", list(ping = TRUE))
  })
}
