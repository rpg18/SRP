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
library(scater)
library(SC3)
library(SingleCellExperiment)
library(plotly)

#load a SingleCellExperiment R object produced by the PigeonPigs pipeline.
#FinalObject includes a pre-run tSNE , PCA and DiffusionMap
load("./data/FinalObject")

#extract required frames from SingleCellExperiment object
TSNEframe<-reducedDim(FinalObject,"TSNE")
PCAframe<-reducedDim(FinalObject,"PCA")

# Define UI for application
ui <- fluidPage(
   
   # Application title
   titlePanel("3D tSNE or PCA plots of PigeonPigs analysis"),
   
   # Sidebar with options 
   sidebarLayout(
      sidebarPanel(
        #select which graph type
        selectInput("Graph1", label = h3("Select graph type"),
                    choice = list("tSNE","PCA"), selected = "PCA"),
        selectInput("Graph2", label = h3("Select graph type"),
                    choice = list("tSNE","PCA"), selected = "PCA"),
        #select colours for graphs
        selectInput("G1colour", label = h3("Select attribute to colour by"), 
                                      choices = list("age", "cell_type", "sc3_2_clusters" , "sc3_3_clusters","sc3_4_clusters","sc3_5_clusters",
                                                     "sc3_6_clusters","sc3_7_clusters","sc3_8_clusters",
                                                     "sc3_9_clusters","sc3_10_clusters","sc3_11_clusters",
                                                     "sc3_12_clusters","sc3_13_clusters","sc3_14_clusters",
                                                     "sc3_15_clusters"),
                                      selected = "cell_type"),
        selectInput("G2colour", label = h3("Select attribute to colour by"), 
                    choices = list("age", "cell_type", "sc3_2_clusters" , "sc3_3_clusters","sc3_4_clusters","sc3_5_clusters",
                                   "sc3_6_clusters","sc3_7_clusters","sc3_8_clusters",
                                   "sc3_9_clusters","sc3_10_clusters","sc3_11_clusters",
                                   "sc3_12_clusters","sc3_13_clusters","sc3_14_clusters",
                                   "sc3_15_clusters"),
                    selected = "cell_type")
      ), 
      # Show plots
      mainPanel(
         plotlyOutput("Grph1"),
         plotlyOutput("Grph2")
      )
   )
)

# Define server logic
server <- function(input, output) {
   #colour_option <- input$PCAcolour
  
   output$Grph1 <- renderPlotly({
     switch(input$Graph1, "tSNE" = plot_ly(as.data.frame(TSNEframe), x = TSNEframe[,1], y= TSNEframe[,2], z= TSNEframe[,3], 
                                           type = 'scatter3d', color = FinalObject[[input$G1colour]], marker = list(size=3)),
     "PCA" = plot_ly(as.data.frame(PCAframe), x = PCAframe[,1], y= PCAframe[,2], z= PCAframe[,3], 
                     type = 'scatter3d', color = FinalObject[[input$G1colour]], marker = list(size=3))
     )
   })
   output$Grph2 <- renderPlotly({
     switch(input$Graph2, "tSNE" = plot_ly(as.data.frame(TSNEframe), x = TSNEframe[,1], y= TSNEframe[,2], z= TSNEframe[,3], 
                                           type = 'scatter3d', color = FinalObject[[input$G2colour]], marker = list(size=3)),
            "PCA" = plot_ly(as.data.frame(PCAframe), x = PCAframe[,1], y= PCAframe[,2], z= PCAframe[,3], 
                            type = 'scatter3d', color = FinalObject[[input$G2colour]], marker = list(size=3))
            )
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

