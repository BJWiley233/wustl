library(shiny)

server <- function(input, output, session) {
  list_values <- c(1:3)
  names(list_values) <- c("One", "Two", "Three")
  rv = reactiveVal(1)
  observeEvent(input$si_page_one_selector,{
    tmp = rv() # Read reactiveVal
    tmp = input$si_page_one_selector
    rv(tmp) # set reactiveVal to new value
  })
  
  observeEvent(input$si_page_two_selector,{
    tmp = rv()
    tmp = input$si_page_two_selector
    rv(tmp)
  })
  # 1. Selector for 'Page one'
  output$uo_page_one_selector <- renderUI({
    # 3. Result
    selectizeInput('si_page_one_selector', "Selector",
                   choices = list_values, multiple = FALSE, 
                   selected = list_values[as.numeric(rv())])
  })
  
  # 2. Selector for 'Page two'
  output$uo_page_two_selector <- renderUI({
    # 3. Result
    selectizeInput('si_page_two_selector', "Selector",
                   choices = list_values, 
                   multiple = FALSE, 
                   selected = list_values[as.numeric(rv())])
  })
  
}

ui <- fluidPage(
  navbarPage("Demo",
             tabPanel("One",
                      fluidPage(
                        fluidRow(uiOutput("uo_page_one_selector")),
                        fluidRow("Page one"))),
             
             tabPanel("Two",
                      fluidPage(
                        fluidRow(uiOutput("uo_page_two_selector")),
                        fluidRow("Page two")))
  )
)

# Run the application 
shinyApp(ui, server)