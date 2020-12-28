list.of.packages <- c("shiny", "ggplot2", "DT", "dplyr", "data.table", "DESeq2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(data.table)
library(DESeq2)



setwd("/home/coyote/Files_for_Oscar/Task09")
## for sourcing the counts and results
##    1. Results for controlling for gender
##    2. Results for controlling for age
##    3. Results for control for both gender and age
#source("task09.R")
#genes <- as.data.table(rownames(count.mat))

## faster to just read in final data that was written than to resource and run DESeq()
## these all come from script 'task09.R'
cleaned_counts <- read.table("counts_cleaned.txt", sep = "\t")
mangled.genes <- read.table("ExcelMangledGenes.txt", skip = 3, header = F)
genes <- as.data.table(row.names(cleaned_counts))
pheno.dat <- read.csv("mayo.path_aging.con.phenotype.csv")
all(sub("^X", "", colnames(cleaned_counts)) == pheno.dat$UID)
res.gender.order.by.fc <- read.table("controlled_gender_res.txt", sep = "\t")
res.age.order.by.fc <- read.table("controlled_age_res.txt", sep = "\t")
res.gender.age.order.by.fc <- read.table("controlled_gender_age_res.txt", sep = "\t")

dseq.gender <- DESeqDataSetFromMatrix(round(cleaned_counts),
                                      pheno.dat,
                                      design = ~ Sex + Diagnosis)
pheno.dat2 <- pheno.dat
pheno.dat2$AgeAtDeath <- as.numeric(gsub("_.*", "", pheno.dat2$AgeAtDeath))
dseq.age <- DESeqDataSetFromMatrix(round(cleaned_counts),
                                   pheno.dat2, 
                                   design = ~ AgeAtDeath + Diagnosis)
dseq.gender.age <- DESeqDataSetFromMatrix(round(cleaned_counts),
                                          pheno.dat2,
                                          design = ~ Sex + AgeAtDeath + Diagnosis)

#####################################################
## UI
#####################################################
ui <- shinyUI(
  fluidPage(theme = "bootstrap.min.css",
            tags$style(HTML("
                            .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter {
                              color: #4d4645;
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
                                                selectizeInput("gene1", "Select a gene for charts:",
                                                               choices = NULL,
                                                               width = "200px")
                                                )
                          ),
                          fluidRow(br()),
                          fluidRow(column(5, 
                                          HTML("<div class='myDiv'>Controlling for <u>Gender</u> (sorted by abs. val. FC)</div>"),                                          
                                          HTML("<br><p style='font-size: 12px'>** The reference is 'Normal' so FC > 0 is positive for 'Pathologic Aging'</p></br>"),
                                          div(DT::dataTableOutput("gender_table"), style = "font-size: 12px;"),
                                          br(),
                                          plotOutput("gender_plot", height = "250px"),
                                          textOutput("test")),
                                   column(1),
                                   column(5, 
                                          HTML("<div class='myDiv'>Controlling for <u>Age</u> (sorted by abs. val. FC)</div>"),
                                          HTML("<br><p style='font-size: 12px'>** The reference is 'Normal' so FC > 0 is positive for 'Pathologic Aging'</p></br>"),
                                          div(DT::dataTableOutput("age_table"), style = "font-size: 12px;"),
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
                                                               width = "200px"))
                          ),
                          fluidRow(column(12, 
                                          HTML("<div class='myDiv'>Controlling for <u>Gender</u> and <u>Age</u> (sorted by abs. val. FC)</div>"),                                          
                                          HTML("<br><p style='font-size: 12px'>** The reference is 'Normal' so FC > 0 is positive for 'Pathologic Aging'</p></br>"),
                                          div(DT::dataTableOutput("gender_age_table"), style = "font-size: 16px;"),
                                          br(),
                                          plotOutput("gender_age_plot", height = "450px")))
                          )),
               tabPanel(span("ExcelMangledGenes", title = "This is the list that converts gene names back from Excel dates"),
                        fluidRow(column(7,HTML(
                          "This is post from <a href='https://www.biostars.org/p/183018/' target='_blank'>Biostars</a> 
                           indicating how storing gene 
                          <br>annotation files in Excel turns genes to Dates."
                        ))),
                        fluidRow(br(), br()),
                        column(6, div(DT::dataTableOutput("mangled"), style = "font-size: 14;"))
                        )
               )
  )
)

#####################################################
## Server
#####################################################
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
  
  
  # https://github.com/rstudio/shiny/issues/1182
  updateSelectizeInput(session = session, inputId = "gene1", choices = genes$V1,
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
  
  ## https://community.rstudio.com/t/how-to-sync-selectizeinput-between-menu-items-for-navbarpage-in-r-shiny/38851
  ## don't think this works
  vals <- reactiveValues(sync = 1)
  observe({
    req(input$gene1)
    vals$sync <- input$gene1
  })
  
  observe({
    req(input$gene2)
    vals$sync <- input$gene2
  })
  
  observe({
    req(vals$sync)
    updateSelectizeInput(session, 'gene1', selected = vals$sync)
    updateSelectizeInput(session, 'gene2', selected = vals$sync)
  })
    
  
  #####################################################
  ## Control for Gender
  #####################################################
  output$gender_table <- DT::renderDataTable({
    
    df <- data.frame(res.gender.order.by.fc)
    df[,1:4] <- signif(df[,1:4], 2)
    df[,5:6] <- signif(df[,5:6], 4)
    if (nchar(input$gene1)) {
      df <- df[input$gene1,]
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
  ## plot for gender
  plot.gender <- reactive({
    req(input$gene1)
  
    p.gender <- plotCounts(dseq.gender, input$gene1, intgroup = c("Diagnosis", "Sex"),
                           returnData = T, normalized = F, transform = F, pc = 1)
    levels(p.gender$Sex) <- c("Female", "Male")
    ggplot(p.gender, aes(x=Diagnosis, y=log2(count), group = 1)) +
      geom_point(position=position_jitter(w=0.1,h=0)) + 
      stat_summary(fun = "mean", geom = "line", color = "red") +
      facet_wrap(~ Sex) +
      ylab("log2(counts + 1)") +
      ggtitle(input$gene1) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size=11),
            strip.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12))
  })
  
  output$gender_plot <- renderPlot({
    plot.gender()
  })
  
  
  #####################################################
  ## Control for Age
  #####################################################
  output$age_table <- DT::renderDataTable({
    
    df <- data.frame(res.age.order.by.fc)
    df[,1:4] <- signif(df[,1:4], 2)
    df[,5:6] <- signif(df[,5:6], 4)
    if (nchar(input$gene1)) {
      df <- df[input$gene1,]
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
  ## plot for age
  plot.age <- reactive({
    req(input$gene1)
    
    p.age <- plotCounts(dseq.age, input$gene1, intgroup = c("Diagnosis", "AgeAtDeath"),
                        returnData = T, normalized = F, transform = F, pc = 1)
    ggplot(p.age) +
      ## ticks are not log2 - https://support.bioconductor.org/p/105938/
      geom_point(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis)) +
      stat_smooth(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis), method = "lm",
                  formula = y ~ x + I(x^2), size = 1, alpha=0.25) +
      ylab("log2(counts + 1)") +
      ggtitle(input$gene1) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size=11),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12))
  })
  
  output$age_plot <- renderPlot({
    plot.age()
  })
  
  
  #####################################################
  ## Control for Gender and Age
  #####################################################
  output$gender_age_table <- DT::renderDataTable({
    
    df <- data.frame(res.gender.age.order.by.fc)
    df[,1:4] <- signif(df[,1:4], 2)
    df[,5:6] <- signif(df[,5:6], 4)
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
  ## plot for gender and age
  plot.gender.age <- reactive({
    req(input$gene2)
    
    p.gender.age <- plotCounts(dseq.gender.age, input$gene2, intgroup = c("Diagnosis", "Sex", "AgeAtDeath"),
                           returnData = T, normalized = F, transform = F, pc = 1)
    levels(p.gender.age$Sex) <- c("Female", "Male")
    ggplot(p.gender.age) +
      ## ticks are not log2 - https://support.bioconductor.org/p/105938/
      geom_point(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis)) +
      stat_smooth(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis), method = "lm",
                  formula = y ~ x + I(x^2), size = 1, alpha=0.25) +
      facet_wrap(~ Sex) +
      ylab("log2(counts + 1)") +
      ggtitle(input$gene2) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size=11),
            strip.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12))
  })
  
  output$gender_age_plot <- renderPlot({
    plot.gender.age()
  })
  
  ## genes to dates and back
  output$mangled <- renderDataTable(
    datatable(mangled.genes, style = "bootstrap", class = "compact",
              colnames = c('Gene', "Excel Date"))
  )

})

shinyApp(ui = ui, server = server)

