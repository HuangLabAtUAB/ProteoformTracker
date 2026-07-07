fluidPage(
  titlePanel("ProteoformTracker"),
  p("Phase 1: core data model and mass engine. UI for scoring/visualization arrives in later phases."),
  h4("Try it: bare-sequence mass lookup"),
  textInput("demo_sequence", "Amino acid sequence", value = "MAGCK"),
  actionButton("demo_calc", "Calculate mass"),
  verbatimTextOutput("demo_result")
)
