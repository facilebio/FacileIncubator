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
  results_stack <- reactiveVal(tibble::tibble(id = character(), analysis = character(), result = list()))
  
  observe({
    shinyjs::disable("remove_module")
  })
  
  ## UI only works with debug = TRUE for some reason (!)
  debug <- TRUE
  bs4dash <- getOption("facile.bs4dash")
  options(facile.bs4dash = FALSE)
  on.exit(options(facile.bs4dash = bs4dash))
  
  x <- eventReactive(input$dataset, {
    switch(req(input$dataset),
           "TCGA" = FacileData:::exampleFacileDataSet(),
           results_stack()[results_stack()$id == input$dataset, "result", drop = TRUE][[1]]
    )
  })
  
  analysisModule <- reactive({
    switch(req(input$analysis),
           "filter" = FacileShine::filteredReactiveFacileDataStore,
           "fdge" = FacileAnalysis::fdgeAnalysis,
           "fpca" = FacileAnalysis::fpcaAnalysis,
           "ffsea" = FacileAnalysis::ffseaAnalysis,
           NULL
    )
  })

  analysisUI <- reactive({
    switch(req(input$analysis),
           "filter" = FacileShine::filteredReactiveFacileDataStoreUI,
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
    
    ui.content <- analysisUI()(ifelse(isolate(req(input$analysis)) == "filter","ds","analysis"), debug = debug)
    
    ui <- tagList(
      ui.content
    )
    
    ui_with_id <- module_UI(module_id, ui)
    
    ## NOTE: immediate = TRUE is necessary!
    insertUI("#gadget_container", "afterEnd", ui_with_id, immediate = TRUE)
    
    shinyjs::disable("add_module")
    shinyjs::enable("remove_module")
    
    isolate(module_res())
  })
  
  observeEvent(input$remove_module, {
    if (NROW(req(module_stack())) > 0) {
      if (debug) {
        print(paste0("removing module ", module_stack()[1]))
      }
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
    req(input$dataset)
    req(input$analysis != "none")
    
    .x <- x()
    fds. <- FacileData::fds(.x)
    samples. <- dplyr::collect(FacileData::samples(.x), n = Inf)
    
    checkmate::assert_class(fds., "FacileDataStore")
    checkmate::assert_class(samples., "facile_frame")

    FacileShine::ReactiveFacileDataStore(fds., "ds", samples = samples.)
  })
  
  module_res <- reactive({
    res <- callModule(analysisModule(), 
               id = ifelse(req(input$analysis) == "filter", "ds", "analysis"),
               rfds = rfds(), 
               aresult = x(), 
               gdb = reactive({sparrow::exampleGeneSetDb()}), 
               path= reactive(rfds()[["parent.dir"]]),
               debug = debug
    )
    if (req(input$analysis) == "filter") {
      return(rfds())
    } else {
      return(res)
    }
  })
  
  observeEvent(input$add_module, {
    input$dataset
    req(input$analysis != "none")
    module_res()
  })
  
  observeEvent(input$remove_module, {
    
    if (is(req(module_res()), "ReactiveFacileDataStore")) {
      res <- module_res()
    } else {
      res <- FacileAnalysis::faro(module_res())
    }
    
    results_stack(
      rbind(
        results_stack(), 
        tibble::tibble(id = paste0("result_", input$add_module), analysis = input$analysis, result = list(res))
      )
    )
    if (debug) {
      print("results stack:")
      print(results_stack())
    }
    
    shinyWidgets::updatePickerInput(session, "dataset", choices = c("TCGA", results_stack()$id), selected = input$dataset)
  })
  
  output$results_list <- renderTable({
    results_stack()[, c("id", "analysis")]
  })
  
})

