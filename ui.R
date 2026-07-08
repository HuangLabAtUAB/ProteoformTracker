fluidPage(
  titlePanel("ProteoformTracker"),
  p("Phase 1: core data model and mass engine. UI for scoring/visualization arrives in later phases."),
  h4("Try it: bare-sequence mass lookup"),
  textInput("demo_sequence", "Amino acid sequence", value = "MAGCK"),
  actionButton("demo_calc", "Calculate mass"),
  verbatimTextOutput("demo_result"),
  hr(),
  h4("Try it: confounding-protein search"),
  p("Reference-proteome mass index: human, reviewed canonical (UniProt UP000005640)."),
  numericInput("confound_mass", "Target intact mass (Da)", value = 45319, min = 0),
  numericInput("confound_window", "Window (± Da)", value = 5, min = 0),
  actionButton("confound_search", "Search"),
  tableOutput("confound_result")
)
