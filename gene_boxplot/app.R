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
library(ggplot2)

#All is a normalised count matrix of only MHCI pathway genes in 466 cells
#All is derived from Darmanis et al. count matrix (see original_data_analysis.R)
load("./data/All")
#MMbox is a normalised count matrix of only MHCI pathway genes in 465 cell
#MMbox is derived from PigeonPigs pipeline analysis of Darmanis et al. raw data
load("./data/MMbox")

# Define UI for application
ui <- fluidPage(
   
   # Application title
   titlePanel("Box plots of expression of major genes of the MHC I pathway"),
   sidebarLayout(
     sidebarPanel(
       selectInput("Gene", label = "Choose Gene", 
              choices = list("HLA-A ", "HLA-B ", "HLA-C ", "B2M ", "TAPBP ", "CALR ", "ERAP1 ",
                             "HSPA5 ", "PDIA3 ", "TAP2 ", "SEC61A1 ", "SEC61A2 ", "SEC61B ", "SEC61G ", "SEC61"))
     ),
   
     
      # Show plots
      mainPanel(
         titlePanel("Boxplot from Darmanis et al. count matrix"),
         plotOutput("Boxplot"),
         titlePanel("Boxplot from PigeonPigs count matrix"),
         plotOutput("Boxplot2")
      )
   )
)

# Define server logic required
server <- function(input, output) {
  
   output$Boxplot <- renderPlot({
     HLA_Box <- ggplot(data=All, aes(x=Sample, y=All[[input$Gene]])) + geom_boxplot()
     abc<-HLA_Box + geom_jitter(shape = 19, position = position_jitter(width = 0.3, height = 0, seed = 2))+ labs(title=input$Gene,x="Type of cell",y="LogCPM")
     plot(abc)
   })
   output$Boxplot2 <- renderPlot({
     widget <- gsub(" ","",input$Gene)
     HLA_Box <- ggplot(data=MMbox, aes(x=Sample, y=MMbox[[widget]])) + geom_boxplot()
     abc<-HLA_Box + geom_jitter(shape = 19, position = position_jitter(width = 0.3, height = 0, seed = 2))+ labs(title=widget,x="Type of cell",y="LogCPM")
     plot(abc)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

