list.of.packages <- c("shiny", "ggplot2", "DT", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(shiny)
library(DT)
library(ggplot2)
library(dplyr)


data("iris")
## K-means unsupervised so don't really need y
iris.y <- iris$Species
iris.X = iris[, -ncol(iris)]

iris.pca <- prcomp(iris.X)
## first 2 PCs
iris.pca2 <- iris.pca$x[,1:2]

ui <- shinyUI(
  fluidPage(theme = "slate.min.css",
            tags$style(HTML("
                            .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter {
                              color: #a9a8ae;
                            }
                            ")),
    titlePanel(title = "Iris K-means"),
    sidebarLayout(
      sidebarPanel(column(5, sliderInput("clusters", "Select K-means clusters",
                            min = 1, max = 5, value = 3)),
                   fluidRow(br()),
                   fluidRow(br()),
                   DT::dataTableOutput("table"), width = 5),
      mainPanel(plotOutput("plot1"), width = 6))
  )
)


server <- shinyServer(function(input, output, session) {
  

  
  observe({
    cls <- kmeans(iris.pca2, centers = input$clusters, iter.max = 20)
    # max allowed 5 clusters so 5 colors
    cls$centers
    cls$cluster
    ## save new dataframe with predicted labels and original attributes
    iris.pca2.labels <- data.frame(iris.pca2, factor(cls$cluster), iris.X)
    ## get average per predicted species across all original attributes
    summary <- iris.pca2.labels %>%
      select(-c("PC1", "PC2")) %>%
      group_by(factor.cls.cluster.) %>%
      summarise(across(everything(), list(mean)))
    summary[,2:5] <- apply(summary[2:5], 2, function(x) round(x, digits=2))
    summary[,1] <- paste("Predicted species", summary$factor.cls.cluster.)

    
    output$plot1 <- renderPlot({
      
      ggplot(data = iris.pca2.labels, aes(PC1, PC2, color = factor.cls.cluster.)) +
        geom_point() +
        scale_colour_discrete(name = "Predicted species") +
        ggtitle(paste("K-means clustering with", input$clusters, "clusters")) +
        theme(plot.title = element_text(hjust = 0.5))
  
    }, bg="transparent")
    
    output$table <- DT::renderDataTable(
      datatable(style = "bootstrap", class = "compact",
                summary,
                selection=list(mode = "single", target = "cell"),
                options = list(
                  initComplete = JS(
                    "function(settings, json) {",
                    "$(this.api().table().header()).css({'color': '#fff'});",
                    "}")),
                colnames = c("Species", "Avg. Sepal Length", "Avg. Sepal Width",
                             "Avg. Petal Length", "Avg. Petal Width")
                )
    )
  })

})

shinyApp(ui = ui, server = server)