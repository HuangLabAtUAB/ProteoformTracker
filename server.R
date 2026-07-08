function(input, output, session) {
  observeEvent(input$demo_calc, {
    output$demo_result <- renderPrint({
      req(input$demo_sequence)
      p <- build_proteoform(
        id = "demo",
        raw_sequence = toupper(input$demo_sequence),
        provenance = "manual"
      )
      m <- proteoform_mass(p)
      list(
        mature_sequence = p$sequence,
        processing_notes = p$metadata$processing_notes,
        monoisotopic_mass = m$mass
      )
    })
  })

  observeEvent(input$confound_search, {
    output$confound_result <- renderTable({
      req(input$confound_mass, input$confound_window)
      if (is.null(reference_mass_index)) {
        return(data.frame(message = "No reference mass index loaded -- run scripts/build_reference_proteome.R"))
      }
      hits <- query_confounding_proteins(
        reference_mass_index,
        target_mass = input$confound_mass,
        window_da = input$confound_window
      )
      if (nrow(hits) == 0) {
        return(data.frame(message = "No entries within this mass window"))
      }
      hits[, c("id", "length", "mass")]
    })
  })
}
