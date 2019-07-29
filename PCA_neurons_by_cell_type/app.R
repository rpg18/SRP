#import libraries.
library(shiny)
library(plotly)

#gre is a PCA object of adult and fetal neuronal cells
#gre was derived from Darmanis et al. count matrix (see original_data_analysis.R)
gre <- readRDS("./data/gre")
#gre is a vector of colours, matched to Darmanis et al. cell type assignment for cells in gre
#collist was derived from Darmanis et al. count matrix (see original_data_analysis.R)
collist <- readRDS("./data/collist")

# Define UI ----
ui <- fluidPage(
titlePanel("PCA plot of adult and fetal neuronal cells"),

  sidebarPanel(
	#Input goes here
    radioButtons("PC1", label = "PCs to plot", choices = list("PC1","PC2","PC3","PC4","PC5"),
                 selected = "PC1"),
    radioButtons("PC2", label = "PCs to plot", choices = list("PC1","PC2","PC3","PC4","PC5"),
                 selected = "PC2"),
    radioButtons("PC3", label = "PCs to plot", choices = list("PC1","PC2","PC3","PC4","PC5"),
                 selected = "PC3")
    ),

  mainPanel(
    plotlyOutput('plot')
  )
)
# Define server logic ----
server <- function(input, output) {
  output$plot <- renderPlotly({
    Xaxis <- switch(input$PC1, "PC1" = 1, "PC2" = 2, "PC3" = 3, "PC4" = 4, "PC5" = 5)
    Yaxis <- switch(input$PC2, "PC1" = 1, "PC2" = 2, "PC3" = 3, "PC4" = 4, "PC5" = 5)
    Zaxis <- switch(input$PC3, "PC1" = 1, "PC2" = 2, "PC3" = 3, "PC4" = 4, "PC5" = 5)
    plot_ly(as.data.frame(gre$ind$coord), 
            x = gre$ind$coord[,Xaxis], y= gre$ind$coord[,Yaxis], z= gre$ind$coord[,Zaxis], 
            type = 'scatter3d', marker= list(color = collist, size=3))
	#red = neurons, green = fetal quiscent, purple = fetal replicating
  })

}

# Run the app ----
shinyApp(ui = ui, server = server)
