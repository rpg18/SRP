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
library(destiny)

#load a SingleCellExperiment R object produced by the PigeonPigs pipeline.
#FinalObject includes a pre-run tSNE , PCA and DiffusionMap
load("./data/FinalObject")

# Define UI for application that draws a graph
ui <- fluidPage(
   
   # Application title
   titlePanel("2D Diffusion Map, PCA or tSNE sized by gene expression"),
   
   # Sidebar with inputs
   sidebarLayout(
      sidebarPanel(
        radioButtons("type", label = h3("Type of dimensional reduction"), 
                     choices = list("PCA","tSNE","Diffusion Map"), 
                    selected = "Diffusion Map"),
        selectInput("cluster", label = h3("Select attribute to colour by"), 
                               choices = list("cell_type","age","sc3_2_clusters" , "sc3_3_clusters","sc3_4_clusters","sc3_5_clusters",
                                              "sc3_6_clusters","sc3_7_clusters","sc3_8_clusters",
                                              "sc3_9_clusters","sc3_10_clusters","sc3_11_clusters",
                                              "sc3_12_clusters","sc3_13_clusters","sc3_14_clusters",
                                              "sc3_15_clusters"),
                               selected = "sc3_9_clusters"),
        textInput("gene", label = h3("Input gene"), value = "PLP1")
      ),
     
      
      # Show a plot
      mainPanel(
         plotOutput("Map")
      )
   )
)

# Define server logic
server <- function(input, output) {
   
   output$Map <- renderPlot({
     widget <- toupper(input$gene)
     switch(input$type, 
            "Diffusion Map" = plotDiffusionMap(FinalObject, colour_by = input$cluster, size_by = widget),
            "PCA" = plotPCA(FinalObject, colour_by = input$cluster, size_by = widget),
            "tSNE" = plotTSNE(FinalObject, colour_by = input$cluster, size_by = widget)
     )
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

