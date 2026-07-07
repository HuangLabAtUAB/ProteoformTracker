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
}
