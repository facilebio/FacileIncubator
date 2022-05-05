library(shiny)

# https://github.com/Appsilon/dynamic-shiny-modules/blob/3b05aad99f633103788b62a94d8ed198ce4b977b/after.R
remove_shiny_inputs <- function(id, .input) {
  invisible(
    lapply(grep(id, names(.input), value = TRUE), function(i) {
      .subset2(.input, "impl")$.values$remove(i)
    })
  )
}

remove_observers <- function(id, .session) {
  invisible(
    lapply(grep(paste0(id, "_observer"), names(.session$userData), value = TRUE),
           function(i) {
             .subset2(.session$userData, i)$destroy()
           })
  )
}

module_UI <- function(id, ui) {
  ns <- NS(id)
  div(id = id, ui)
}

shinyServer(function(input, output, session) {
  
  module_stack <- reactiveVal(NULL)
  results_stack <- reactiveVal(data.frame(id = integer(), analysis = character()))
  
  observe({
    shinyjs::disable("remove_module")
  })
  
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
    
  observeEvent(input$add_module, {
    req(input$analysis != "none")
    
    # store the id of the newly added module using the 
    # value of the actionButton to make it unique
    module_id <- paste0("id_", input$add_module)
    if (debug) print(paste0("this module is ", module_id))
    module_stack(c(module_id, module_stack()))
    
    ui.content <- analysisUI()("analysis", debug = debug)
    
    ui <- tagList(
      FacileShine::filteredReactiveFacileDataStoreUI("ds"),
      tags$hr(),
      ui.content
    )
    
    ui_with_id <- module_UI(module_id, ui)
    
    ## NOTE: immediate = TRUE is necessary!
    insertUI("#gadget_container", "afterEnd", ui_with_id, immediate = TRUE)
    
    shinyjs::disable("add_module")
    shinyjs::enable("remove_module")
  })
  
  observeEvent(input$remove_module, {
    if (length(module_stack()) > 0) {
      if (debug) print(paste0("removing module ", module_stack()[1]))
      removeUI(paste0("#", module_stack()[1]))
    }
    remove_shiny_inputs(module_stack()[1], input)
    remove_observers(module_stack()[1], session)
    module_stack(module_stack()[-1])
    
    shinyjs::enable("add_module")
    shinyjs::disable("remove_module")
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
    
    FacileShine::ReactiveFacileDataStore(fds., "ds", samples = samples.)
  })
  
  observeEvent(req(input$add_module), {
    req(input$analysis != "none")
    if (debug) print("running analysis")
    analysis <- callModule(analysisModule(), "analysis", rfds(), debug = debug)
    results_stack(rbind(results_stack(), data.frame(id = input$add_module, type = input$analysis)))
  })
  
  output$results_list <- renderTable({
    if (NROW(results_stack()) > 0) results_stack() else NULL
  })
  
})

