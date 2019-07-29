#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#import libraries.
library(shiny)
library(SingleCellExperiment)
library(SC3)
library(scater)

#load a SingleCellExperiment R object produced by the PigeonPigs pipeline.
#FinalObject includes a pre-run tSNE , PCA and DiffusionMap
load("./data/FinalObject")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("PigeonPigs clusters fitting to age and cell type"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
     sidebarPanel(selectInput("cluster", label = h3("Select attribute to colour by"), 
                              choices = list("sc3_2_clusters" , "sc3_3_clusters","sc3_4_clusters","sc3_5_clusters",
                                             "sc3_6_clusters","sc3_7_clusters","sc3_8_clusters",
                                             "sc3_9_clusters","sc3_10_clusters","sc3_11_clusters",
                                             "sc3_12_clusters","sc3_13_clusters","sc3_14_clusters",
                                             "sc3_15_clusters"),
                              selected = "sc3_9_clusters")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("graph")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$graph <- renderPlot({
     plotScater(FinalObject, block1 = "age", block2 = "cell_type",
                colour_by = input$cluster, nfeatures = 200, exprs_values = "counts")
   }, height = 800)
}

# Run the application 
shinyApp(ui = ui, server = server)

