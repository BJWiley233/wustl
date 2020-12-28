list.of.packages <- c("shiny", "ggplot2", "DT", "dplyr", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(data.table)



setwd("/home/coyote/Files_for_Oscar/Task09")
## for sourcing the counts and results
##    1. Results for controlling for gender
##    2. Results for controlling for age
##    3. Results for control for both gender and age
#source("task09.R")
genes <- as.data.table(rownames(count.mat))


ui <- shinyUI(
  fluidPage(theme = "bootstrap.min.css",
            tags$style(HTML("
                            .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter {
                              color: #4d4645;
                            }
                            
                            .dataTables_wrapper {
                              font-size: 12px;
                            }
                            
                            body {
                              background: linear-gradient(120deg, #f75724, #ECE3E6);
                              width: 100%;
                              height: 100%;
                            }

                            .selectize-input > input[placeholder] {
                              width: 100% !important;
                            }
                            
                            a {
                              color: blue
                            }
                            
                            .myDiv {
                              text-align: center;
                               font-weight: bold;
                               font-size: 16px;
                            }
                            
                            .pagination > li > a, .pagination > li > span {
                              padding: 4px 6px;
                            }
                            
                            ")),
    fluidRow(br()),
    navbarPage("Transcriptomic Analysis",
               navbarMenu("Differential Expression",
                 tabPanel(title = "Controlling for Gender & Age Individually",
                          fluidRow(column(7,HTML(
                            "This page will report the differentially expressed genes for the
                            <br>
                            condition 'Diagnosis' (Pathological aging vs. Controls) controlling for 
                            <br>
                            gender and age <a href='https://support.bioconductor.org/p/110873/' target='_blank'>individually</a>. To report for controlling for both gender
                            <br>
                            and age <a href='https://support.bioconductor.org/p/63766/' target='_blank'>together</a>, select under 'Controlling for both Gender & Age' 
                            <br>
                            under 'Differential Expression' above.
                            <br>
                            "
                          ))),
                          fluidRow(br(), column(12, align = "center",
                                                selectizeInput("gene", "Select a gene for charts:",
                                                               choices = NULL,
                                                               width = "200px")
                                                )
                          ),
                          fluidRow(br()),
                          fluidRow(column(5, 
                                          HTML("<div class='myDiv'>Controlling for <u>Gender</u> (sorted by abs. val. FC)</div>"),                                          
                                          HTML("<br><p style='font-size: 12px'>** The reference is 'Normal' so FC > 0 is positive for 'Pathologic Aging'</p></br>"),
                                          div(DT::dataTableOutput("gender_table")),
                                          br(),
                                          plotOutput("gender_plot", height = "250px"),
                                          textOutput("test")),
                                   column(1),
                                   column(5, 
                                          HTML("<div class='myDiv'>Controlling for <u>Age</u> (sorted by abs. val. FC)</div>"),
                                          HTML("<br><p style='font-size: 12px'>** The reference is 'Normal' so FC > 0 is positive for 'Pathologic Aging'</p></br>"),
                                          div(DT::dataTableOutput("age_table")),
                                          br(),
                                          plotOutput("age_plot", height = "250px")
                                          ))
                 ),
                 
                 
                 tabPanel(title = "Controlling for both Gender & Age",
                          fluidRow(column(7,HTML(
                            "This page will report the differentially expressed genes for the
                            <br>
                            condition 'Diagnosis' (Pathological aging vs. Controls) controlling for 
                            <br>
                            for both gender and age <a href='https://support.bioconductor.org/p/63766/' target='_blank'>together</a>.
                            <br>
                            "
                          ))),
                          fluidRow(br(), column(12, align = "center",
                                                selectizeInput("gene2", "Select a gene for charts:",
                                                               choices = NULL,
                                                               width = "200px")
                                                )
                          ),
                          
                          
                          fluidRow(column(12, 
                                          HTML("<div class='myDiv'>Controlling for <u>Gender</u> and <u>Age</u> (sorted by abs. val. FC)</div>"),                                          
                                          HTML("<br><p style='font-size: 12px'>** The reference is 'Normal' so FC > 0 is positive for 'Pathologic Aging'</p></br>"),
                                          div(DT::dataTableOutput("gender_age_table")),
                                          br(),
                                          plotOutput("gender_age_plot", height = "450px")))
                          
                          
                          )),
               tabPanel(span("ExcelMangledGenes", title = "This is the list that converts gene names back from Excel dates"),
                        column(12, DT::dataTableOutput("table2"))
                        )
               )
  )
)


server <- shinyServer(function(input, output, session) {
  
  ## https://community.rstudio.com/t/use-of-regular-expressions-in-selectizeinput-in-shiny-app/13064/2
  getScore <- function() {
    return(I("function(search) {
                var score = this.getScoreFunction(search);
                return function(item) {
                  return item.label.toLowerCase().startsWith(search.toLowerCase()) ? 1 : 0;
                };
              }"
      )
    )
  }
  
  rv <- reactiveVal(1)
  
  observeEvent(input$gene, {
    tmp <- rv()
    tmp <- input$gene
    rv(tmp)
    cat(rv)
  })
  
  observeEvent(input$gene2, {
    tmp <- rv()
    tmp <- input$gene2
    rv(tmp)
  })
  
  ## https://github.com/rstudio/shiny/issues/1182
  updateSelectizeInput(session = session, inputId = "gene", choices = genes$V1,
                       server = T, selected = character(0),
                       options = list(
                         placeholder = "Type to select a gene",
                         items = c(),
                         score = getScore(),
                         dropdownParent = 'body',
                         openOnFocus = FALSE
                       ))
  
  updateSelectizeInput(session = session, inputId = "gene2", choices = genes$V1,
                       server = T, selected = character(0),
                       options = list(
                         placeholder = "Type to select a gene",
                         items = c(),
                         score = getScore(),
                         dropdownParent = 'body',
                         openOnFocus = FALSE
                       ))
  
  output$test <- renderText({ !nchar(input$gene) })
  
  ## update each page for the gene selected if you want to change between
  ## controlling for gender/age individually and together
  ## https://stackoverflow.com/questions/57738924/how-to-sync-selectizeinput-between-menu-items-for-navbarpage-in-r-shiny
  # rv <- reactiveVal(1)
  # 
  # observeEvent(input$gene, {
  #  tmp <- rv()
  #  tmp <- input$gene
  #  rv(tmp)
  # })
  # 
  # observeEvent(input$gene2, {
  #   tmp <- rv()
  #   tmp <- input$gene2
  #   rv(tmp)
  # })
  
  #####################################################
  ## Control for Gender
  #####################################################
  output$gender_table <- DT::renderDataTable({
    
    df <- round(data.frame(res.sex.order.by.fc),2)
    if (nchar(input$gene)) {
      df <- df[input$gene,]
    }
    datatable(style = "bootstrap", class = "compact",
              ## https://stackoverflow.com/questions/31124122/r-shiny-mouseover-text-for-table-columns
              callback = JS("
                var tips = ['Gene Name', 'Base Mean', 'Log2 Fold Change, the reference is `Control` so > 0 means + for `Pathologic Aging`', 
                            'Log2 Fold Change Standard Error','Wald statistic', 'p-value', 'FDR adjusted p-value'],
                header = table.columns().header();
                for (var i = 0; i < tips.length; i++) {
                  $(header[i]).attr('title', tips[i]);
                }
                "
              ),
              df, selection=list(mode = "single", target = "row")
    )
    
  })
  
  
  plot.gender <- reactive({
    req(input$gene)
  
    p.gender <- plotCounts(dseq.sex, input$gene, intgroup = c("Diagnosis", "Sex"),
                           returnData = T, normalized = F, transform = F, pc = 1)
    levels(p.gender$Sex) <- c("Female", "Male")
    ggplot(p.gender, aes(x=Diagnosis, y=log2(count), group = 1)) +
      geom_point(position=position_jitter(w=0.1,h=0)) + 
      stat_summary(fun = "mean", geom = "line", color = "red") +
      facet_wrap(~ Sex) +
      ylab("log2(counts + 1)") +
      ggtitle(input$gene) +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$gender_plot <- renderPlot({
    plot.gender()
  })
  
  
  #####################################################
  ## Control for Age
  #####################################################
  output$age_table <- DT::renderDataTable({
    
    df <- round(data.frame(res.age.order.by.fc),2)
    if (nchar(input$gene)) {
      df <- df[input$gene,]
    }
    datatable(style = "bootstrap", class = "compact",
              callback = JS("
                var tips = ['Gene Name', 'Base Mean', 'Log2 Fold Change, the reference is `Control` so > 0 means + for `Pathologic Aging`', 
                            'Log2 Fold Change Standard Error','Wald statistic', 'p-value', 'FDR adjusted p-value'],
                header = table.columns().header();
                for (var i = 0; i < tips.length; i++) {
                  $(header[i]).attr('title', tips[i]);
                }
                "
              ),
              df, selection=list(mode = "single", target = "row")
    )
    
  })
  
  plot.age <- reactive({
    req(input$gene)
    
    p.age <- plotCounts(dseq.age, input$gene, intgroup = c("Diagnosis", "AgeAtDeath"),
                        returnData = T, normalized = F, transform = F, pc = 1)
    ggplot(p.age) +
      ## ticks are not log2 - https://support.bioconductor.org/p/105938/
      geom_point(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis)) +
      stat_smooth(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis), method = "lm",
                  formula = y ~ x + I(x^2), size = 1, alpha=0.25) +
      ylab("log2(counts + 1)") +
      ggtitle(input$gene) +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$age_plot <- renderPlot({
    plot.age()
  })
  
  
  #####################################################
  ## Control for Gender and Age
  #####################################################
  output$gender_age_table <- DT::renderDataTable({
    
    df <- round(data.frame(res.gender.age.order.by.fc),2)
    if (nchar(input$gene2)) {
      df <- df[input$gene2,]
    }
    datatable(style = "bootstrap", class = "compact",
              ## https://stackoverflow.com/questions/31124122/r-shiny-mouseover-text-for-table-columns
              callback = JS("
                var tips = ['Gene Name', 'Base Mean', 'Log2 Fold Change, the reference is `Control` so > 0 means + for `Pathologic Aging`', 
                            'Log2 Fold Change Standard Error','Wald statistic', 'p-value', 'FDR adjusted p-value'],
                header = table.columns().header();
                for (var i = 0; i < tips.length; i++) {
                  $(header[i]).attr('title', tips[i]);
                }
                "
              ),
              df, selection=list(mode = "single", target = "row")
    )
    
  })
  
  
  plot.gender.age <- reactive({
    req(input$gene2)
    
    dseq.gender.age <- plotCounts(dseq.sex, input$gene2, intgroup = c("Diagnosis", "Sex", "AgeAtDeath"),
                           returnData = T, normalized = F, transform = F, pc = 1)
    levels(dseq.gender.age$Sex) <- c("Female", "Male")
    ggplot(p.gender.age) +
      ## ticks are not log2 - https://support.bioconductor.org/p/105938/
      geom_point(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis)) +
      stat_smooth(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis), method = "lm",
                  formula = y ~ x + I(x^2), size = 1, alpha=0.25) +
      facet_wrap(~ Sex) +
      ylab("log2(counts + 1)") +
      ggtitle(input$gene2) +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$gender_age_plot <- renderPlot({
    plot.gender.age()
  })

})

shinyApp(ui = ui, server = server)

