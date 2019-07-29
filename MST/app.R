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
library(igraph)

#fetaltree is an igraph object of fetal cells only
#fetaltree was derived from Darmanis et al. count matrix (see original_data_analysis.R)
load("./data/fetaltree")
#neurontree is an igraph object of adult neuronal cells only
#neurontree was derived from Darmanis et al. count matrix (see original_data_analysis.R)
load("./data/neurontree")
# Define UI for application
ui <- fluidPage(
  
      # Show plots
      mainPanel(
         titlePanel("Minimum Spanning Tree of Fetal Neurons"),
         plotOutput("MST"),
         titlePanel("Minimum Spanning Tree of Adult Neurons"),
         plotOutput("MST2")
      )
)

# Define server logic
server <- function(input, output) {
   
   output$MST <- renderPlot({
     plot(fetreegra, layout=layout.fruchterman.reingold,
          vertex.size=4, vertex.label=NA, asp=FALSE, edge.arrow.mode=0)
   })
   output$MST2 <- renderPlot({
     plot(nertreegra, layout=layout.fruchterman.reingold,
          vertex.size=4, vertex.label=NA, asp=FALSE, edge.arrow.mode=0)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

