library(shiny)

shinyUI(fluidPage(

    # Application title
    titlePanel("FacileSandbox"),

    sidebarLayout(
        sidebarPanel(width = 3,
                     shinyjs::useShinyjs(),
                     shinyWidgets::pickerInput("dataset",
                                               "Select an input:",
                                               choices = c("TCGA"),
                                               selected = "TCGA"
                     ),
                     shinyWidgets::pickerInput("analysis",
                                               "Select an output:",
                                               choices = c("none", "filter", "fdge", "fpca", "ffsea"),
                                               selected = "none"
                     ),
                     actionButton("add_module", "Add", icon = icon("plus-circle")),
                     actionButton("remove_module", "Done", icon = icon("check-circle")),
                     br(),
                     h4("Results:"),
                     tableOutput("results_list")
        ),
        
        mainPanel(
          tags$div(id = "gadget_container")
        )
    )
))
