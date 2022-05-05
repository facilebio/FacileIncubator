library(shiny)

shinyUI(fluidPage(

    # Application title
    titlePanel("FacileSandbox"),

    sidebarLayout(
        sidebarPanel(width = 3,
          
          shinyWidgets::pickerInput("dataset",
                                    "Select a dataset",
                                    choices = "TCGA",
                                    selected = "TCGA"
          ),
          shinyWidgets::pickerInput("analysis",
                                    "Select an analysis:",
                                    choices = c("none", "fdge", "fpca", "ffsea"),
                                    selected = "none"
          )
        ),

        mainPanel(
            tags$div(id = "here")
        )
    )
))
