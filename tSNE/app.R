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
library(plotly)

#sne is a table output from R containing the dimensions of a tSNE run and cell ref
#sne is whitespace separated
#sne was derived from Darmanis et al. count matrix (see original_data_analysis.R)
sne <- read.table("./data/sne")
#cluslist is a vector of colours, matched to our cluster assignment for cells in sne
#cluslist was derived from Darmanis et al. count matrix (see original_data_analysis.R)
load("./data/cluslist")
#grahamlist is a vector of colours, matched to Darmanis et al. cell type assignment for cells in sne
#grahamlist was derived from Darmanis et al. count matrix (see original_data_analysis.R)
load("./data/grahamlist")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("tSNE plot"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(radioButtons("colour", label = "Colour by...", choices = list("Cluster","Cell type"),
                                selected = "Cluster")
      ),
      # Show a plot of the generated distribution
      mainPanel(
         plotlyOutput("tSNE")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
   output$tSNE <- renderPlotly({
     widget <- switch(input$colour, "Cluster" = cluslist, "Cell type" = grahamlist)
     plot_ly(as.data.frame(sne), x = sne[,2], y= sne[,3], z= sne[,4], 
             type = 'scatter3d', marker= list(color = widget, size=1)
             )
   })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

