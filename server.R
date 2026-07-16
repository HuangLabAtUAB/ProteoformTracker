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
  # suspendWhenHidden = FALSE treatment.
  output$digestion_candidates_ui <- renderUI(NULL)
  output$digestion_coverage_text <- renderText("")
  for (nm in c("gene_status", "isoform_catalog_ui", "ptm_warnings_ui", "proteoform_table",
               "stale_notice", "viz_script", "confounder_status",
               "digestion_candidates_ui", "digestion_coverage_text")) {
    outputOptions(output, nm, suspendWhenHidden = FALSE)
  }

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
      return()
    }
    catalog(cat_result)

    withProgress(message = paste("Fetching real protein sequences for", gene), value = 0, {
      n <- length(cat_result$transcript_ids)
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
      bare_pfs(pfs)
    })

    n_ok <- sum(!vapply(bare_pfs(), is.null, logical(1)))
    output$gene_status <- renderText(sprintf(
      "%s: %d transcripts in the precomputed exon index, %d with a real fetched protein sequence.",
      gene, length(cat_result$transcript_ids), n_ok
    ))
  })

  output$isoform_catalog_ui <- renderUI({
    cat_result <- catalog()
    pfs <- bare_pfs()
    if (is.null(cat_result) || is.null(pfs)) return(tags$p(class = "pt-note", "Load a gene to see its real isoform catalog."))

    rows <- lapply(cat_result$transcript_ids, function(tid) {
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
      tags$div(class = "pt-isorow",
        checkboxInput(paste0("iso_chk_", tid), NULL, value = default_checked),
        tags$span(class = "pt-id", tid),
        tags$span(class = "pt-meta", sprintf("%d aa, %.1f Da", nchar(pf$sequence), mass)),
        tags$span(class = "pt-meta", sprintf("exons %s", compress_exon_ranges(exon_nums))),
        textInput(paste0("iso_ptm_", tid), NULL, value = ptm_default,
                  placeholder = "e.g. 133_Thr_Phospho; 210_Pro_Oxidation,215_Ser_Sulfo", width = "420px")
      )
    })
    tags$div(rows)
  })

  # Recomputes on ANY iso_chk_*/iso_ptm_* change -- Shiny tracks dynamic
  # input[[...]] reads at runtime, so this reactive correctly re-fires even
  # though the input IDs themselves are generated dynamically above.
  derived <- reactive({
    cat_result <- catalog()
    pfs <- bare_pfs()
    req(cat_result, pfs)
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
    req(length(d$rows) > 0)
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

    header <- tags$div(class = "pt-isorow", style = "font-weight:600;color:#555;font-size:12px;",
      tags$span(style = "width:20px;"), tags$span(style = "min-width:150px;", "Parent"),
      tags$span(style = "min-width:110px;", "Range"), tags$span(style = "min-width:90px;", "Missed cl."),
      tags$span(style = "min-width:100px;", "Mass (Da)"), tags$span(style = "min-width:110px;", "MS1 FWHM (Da)"),
      tags$span(style = "min-width:110px;", "MS2 propensity"), tags$span(style = "min-width:90px;", "PTM sites")
    )
    rows <- lapply(seq_len(nrow(cands)), function(r) {
      cid <- cands$id[r]
      prev_checked <- isolate(input[[paste0("pep_chk_", sanitize_html_id(cid))]])
      default_checked <- if (!is.null(prev_checked)) prev_checked else (cid %in% default_ids)
      tags$div(class = "pt-isorow",
        checkboxInput(paste0("pep_chk_", sanitize_html_id(cid)), NULL, value = default_checked),
        tags$span(class = "pt-meta", style = "min-width:150px;", dg$parent_label[[cid]] %||% cands$parent_id[r]),
        tags$span(class = "pt-meta", style = "min-width:110px;", sprintf("%d-%d (%d aa)", cands$start[r], cands$end[r], cands$length[r])),
        tags$span(class = "pt-meta", style = "min-width:90px;", cands$missed_cleavages[r]),
        tags$span(class = "pt-meta", style = "min-width:100px;", sprintf("%.1f", cands$mass[r])),
        tags$span(class = "pt-meta", style = "min-width:110px;", sprintf("%.2f", cands$ms1_fwhm_da[r])),
        tags$span(class = "pt-meta", style = "min-width:110px;", sprintf("%.2f", cands$ms2_avg_propensity[r])),
        tags$span(class = "pt-meta", style = "min-width:90px;", cands$ptm_sites_covered[r])
      )
    })
    tags$div(header, tags$div(rows))
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
      reference_mass_index
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
    req(catalog())
    "Click \"Run analysis\" to (re)compute MS1/MS2/confounder results for the currently included proteoforms."
  })

  output$viz_script <- renderUI({
    a <- analysis()
    j1 <- jsonlite::toJSON(a$payload1, auto_unbox = TRUE, digits = 4, null = "null")
    j2 <- if (!is.null(a$payload2)) jsonlite::toJSON(a$payload2, auto_unbox = TRUE, digits = 4, null = "null") else "null"
    tags$script(HTML(sprintf(
      "PT.renderSection1(%s); var __pt_p2 = %s; if (__pt_p2) { PT.renderSection2(__pt_p2); } else { document.getElementById('s2-ladder').innerHTML=''; document.getElementById('s2-ms1').innerHTML=''; document.getElementById('s2-legend').innerHTML='<p style=\"color:#888;font-size:12px;\">No confounder-search target selected, or no reference mass index loaded.</p>'; }",
      j1, j2
    )))
  })

  output$confounder_status <- renderText({
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

  observeEvent(input$btn_run_fasta, {
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
    req(a, is.null(a$align_error))

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
      selected_pf <- build_proteoforms_for_transcripts(selected_tids)
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
      catalog(list(gene = gene_name, exon_table = novel_exon_table,
                    transcript_ids = "NOVEL_1", protein_lengths = c(NOVEL_1 = nchar(orf_row$protein_sequence))))
      bare_pfs(setNames(list(novel_pf), "NOVEL_1"))
    }
    precheck_extra(extra_checked)

    output$gene_status <- renderText(sprintf(
      "From FASTA input: %s translation (%d aa) added%s. Pick it (and any known isoforms) below, then Run analysis.",
      orf_row$orf_id, nchar(orf_row$protein_sequence),
      if (!is.na(gene_name)) paste0(" alongside ", gene_name, "'s known isoforms") else ""
    ))
    updateRadioButtons(session, "input_mode", selected = "gene")
  })

  # ============================================================
  # Option 3: rMATS upload -- not yet wired (adapter from IsoPepTracker's
  # rMATS logic is a separate follow-up pass).
  # ============================================================
  observeEvent(input$btn_run_rmats, {
    output$rmats_status <- renderText("rMATS -> exon-structure adapter not yet implemented in this app (planned follow-up).")
  })
}
