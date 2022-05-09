library(shiny)

shinyUI(fluidPage(

    # Application title
    titlePanel("FacileSandbox"),

    sidebarLayout(
        sidebarPanel(width = 3,
                     shinyjs::useShinyjs(),
                     shinyWidgets::pickerInput("dataset",
                                               "Select an input:",
                                               choices = "TCGA",
                                               selected = "TCGA"
                     ),
                     shinyWidgets::pickerInput("analysis",
                                               "Select an output:",
                                               choices = c("none", "fdge", "fpca", "ffsea"),
                                               selected = "none"
                     ),
                     actionButton("add_module", "Add", icon = icon("plus-circle")),
                     actionButton("remove_module", "Remove", icon = icon("trash-alt")),
                     br(),
                     h4("Results:"),
                     tableOutput("results_list")
        ),
        
        mainPanel(
          tags$div(id = "gadget_container")
        )
    )
))
