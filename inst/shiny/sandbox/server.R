library(shiny)

shinyServer(function(input, output, session) {
  
  debug <- FALSE
  bs4dash <- getOption("facile.bs4dash")
  options(facile.bs4dash = FALSE)
  on.exit(options(facile.bs4dash = bs4dash))
  
  x <- reactive({
    switch(req(input$dataset),
           "TCGA" = FacileData:::exampleFacileDataSet(),
           NULL
    )
  })
  
  analysisModule <- reactive({
    switch(req(input$analysis),
           "fdge" = FacileAnalysis::fdgeAnalysis,
           "fpca" = FacileAnalysis::fpcaAnalysis,
           "ffsea" = FacileAnalysis::ffseaAnalysis,
           NULL
    )
  })

  analysisUI <- reactive({
    switch(req(input$analysis),
           "fdge" = FacileAnalysis::fdgeAnalysisUI,
           "fpca" = FacileAnalysis::fpcaAnalysisUI,
           "ffsea" = FacileAnalysis::ffseaAnalysisUI,
           NULL
    )
  })
    
  observeEvent(input$analysis, {
    req(input$analysis != "none")
    ui.content <- analysisUI()("analysis", debug = debug)
    
    ui <- tagList(
      FacileShine::filteredReactiveFacileDataStoreUI("ds"),
      tags$hr(),
      ui.content
    )
    
    ## NOTE: immediate = TRUE is necessary!
    insertUI("#here", "afterEnd", ui, immediate = TRUE)
  })
  
  ## this logic should be isolated into a function
  rfds <- reactive({
    req(input$analysis != "none")
    
    .x <- x()
    if (is(.x, "facile_frame")) {
      fds. <- FacileData::fds(.x)
      samples. <- .x
      sample.filter <- FALSE
      restrict_samples <- samples.
    } else if (is(.x, "FacileDataStore")) {
      sample.filter <- TRUE
      fds. <- .x
      samples. <- dplyr::collect(FacileData::samples(.x), n = Inf)
    } else if (is(.x, "FacileAnalysisResult")) {
      # ugh, this isn't going to work -- I'm writing this in to fire up a
      # ffseaGadget, whose needs to be a FacileAnalysisResult.
      sample.filter <- FALSE
      fds. <- FacileData::fds(.x)
      samples. <- dplyr::collect(FacileData::samples(.x), n = Inf)
      restrict_samples <- samples.
    } else {
      stop("What in the world?")
    }
    checkmate::assert_class(fds., "FacileDataStore")
    checkmate::assert_class(samples., "facile_frame")
    
    FacileShine:::ReactiveFacileDataStore(fds., "ds", samples = samples.)
  })
  
  observe({
    req(input$analysis != "none")
    analysis <- callModule(analysisModule(), "analysis", rfds(), debug = debug)
  })
  
})

