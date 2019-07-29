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
library(mclust)

#load a mclust object derived from a tsne object
#yoghurt is derived from Darmanis et al. count matrix (see original_data_analysis.R)
load("./data/yoghurt")

# Define UI for application
ui <- fluidPage(
   
      # Show plots
      mainPanel(
        titlePanel("Bayesian Information Criterion plot"),
         plotOutput("BIC"),
        titlePanel("Boxplot of uncertainty of classification within members of each cluster"),
         plotOutput("Boxplot")
      )
   )


# Define server logic required
server <- function(input, output) {
   
   output$BIC <- renderPlot({
     plot(yoghurt$BIC, colors = c("gray", "black", "orange", "darkred", "red", "magenta", "darkgreen", "green","lightblue","darkblue")) 
   })
   output$Boxplot <- renderPlot({
     boxplot(yoghurt$z)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

