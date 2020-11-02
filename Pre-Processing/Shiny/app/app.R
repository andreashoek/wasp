## R libraries

library(shinydashboard)
library(shiny)
library(shinyFiles)
library(rjson)
library(DT)
library(ggplot2)
library(readr)
library(plotly)
library(scales)
library(shinyjs)

options(shiny.maxRequestSize=500*1024^2)

## UI
ui <- 
  dashboardPage(
    # 
    dashboardHeader(tags$li(class = "dropdown", fluidRow( shinyjs::useShinyjs(),
                                                         column(7, actionButton("help", label = "Information", icon = icon("fas fa-info-circle")))), position = "fixed-top")),
    dashboardSidebar(
      sidebarMenuOutput("menu")
    ),
    dashboardBody(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
      ),
      tabItems(
        tabItem(tabName = "about",
                fluidRow(
                  column(12,
                         h2("Welcome to WASP!"),
                         br(),
                         p("This is a quality control tool for visualizing Pre-Processing of ddSeq-based single-cell RNA sequencing data.
                                    It is designed to interactively visualize results of the WASP Snakemake pipeline and to help users get insight into their data and validate cells for further downstream analysis."),
                         br(),
                         p("For the proper visualization of the data it is important that the preprocessing was accomplished with the WASP Pre-Processing snakemake pipeline
                                    and the results were saved in JSON files."),
                         hr(style="border-color:black;"),
                         br(),
                         p("Please select the path to your Results directory."),
                         br(),
                         tags$div(tags$label("Choose folder", class="btn btn-primary",
                                             tags$input(id = "fileIn", type = "file", style="display: none;", onchange="pressed()", webkitdirectory = "TRUE"))),
                         br(),
                         tags$div(id="fileIn_progress", class="progress progress-striped active shiny-file-input-progress",
                                  tags$div(class="progress-bar")),
                         br(),
                         h5("Samples"),
                         verbatimTextOutput("path", placeholder = TRUE),
                         h5("Info"),
                         verbatimTextOutput("info"),
                         br(),
                         fluidRow(align="center", uiOutput("welcomeButton"))
                  )
                )
        ),
        tabItem(tabName = "sampleSummaries",
                htmlOutput("sampTop"),
                fluidRow(
                  column(6, align="left", h3("Summary of all samples")),
                  column(6, align="right", actionButton(inputId = "sumButton", label = "Continue"))),
                br(),
                fluidRow(
                  column(12,
                         dataTableOutput("summaryTab"),
                         br(),
                         h3("Counted reads of samples"),
                         plotOutput("sumBarChart"),
                         fluidRow(column(12, align= "right", downloadButton("downloadSampSum", "Download plot"))),
                         br(),
                         h3("Mapping metrics of STAR analysis"),
                         dataTableOutput("mappingSum"),
                         br(),
                         h3("Mapping metrics of samples"),
                         plotOutput("mapSumChart"),
                         fluidRow(column(12, align= "right", downloadButton("downloadSampMap", "Download plot"))),
                         br(),
                         h3("Gene count metrics of featureCount analysis"),
                         dataTableOutput("geneSum"),
                         br(),
                         h3("Gene count metrics of samples"),
                         plotOutput("geneSumChart"),
                         fluidRow(column(12, align= "right", downloadButton("downloadSampGene", "Download plot"))),
                         fluidRow(align="center", htmlOutput("topSum"))
                  )
                )
        ),
        tabItem(tabName = "fastqcResults",
                fluidRow(
                  uiOutput("fcSamples"),
                  column(12,
                         sidebarLayout(
                           sidebarPanel(
                             style = "position:fixed;",
                             uiOutput("basic"),
                             uiOutput("baseSeq"),
                             uiOutput("tileQua"),
                             uiOutput("seqQua"),
                             uiOutput("seqCon"),
                             uiOutput("gcCon"),
                             uiOutput("nCon"),
                             uiOutput("lenDist"),
                             uiOutput("dupLev"),
                             uiOutput("overrepSeq"),
                             uiOutput("adapCon"),
                             width = 3
                           ),
                           mainPanel(
                             fluidRow(
                               column(11, uiOutput("basicSta", style="padding-top: 30px; margin-top: -5px;"),
                                      dataTableOutput("basicstats")),
                               column(1, align="right", actionButton(inputId = "fastButton", label = "Continue"), style= "margin-top:20px;")),
                             br(),
                             uiOutput("perBaseQua", style="padding-top: 40px; margin-top: -40px;"),
                             imageOutput("perBaseQuality", height = "600px"),
                             br(),
                             uiOutput("tileQuali", style="padding-top: 40px; margin-top: -40px;"),
                             imageOutput("tileQuality", height = "600px"),
                             br(),
                             uiOutput("qualitySco", style="padding-top: 40px; margin-top: -40px;"),
                             imageOutput("qualityScores", height = "600px"),
                             br(),
                             uiOutput("sequenceCont", style="padding-top: 40px; margin-top: -40px;"),
                             imageOutput("sequenceContent", height = "600px"),
                             br(),
                             uiOutput("gcCont", style="padding-top: 40px; margin-top: -40px;"),
                             imageOutput("gcContent", height = "600px"),
                             br(),
                             uiOutput("nCont", style="padding-top: 40px; margin-top: -40px;"),
                             imageOutput("nContent", height = "600px"),
                             br(),
                             uiOutput("lengthDist", style="padding-top: 40px; margin-top: -40px;"),
                             imageOutput("lengthDistribution", height = "600px"),
                             br(),
                             uiOutput("duplicationLev", style="padding-top: 40px; margin-top: -40px;"),
                             imageOutput("duplication", height = "600px"),
                             br(),
                             uiOutput("overrepresentedSeq", style="padding-top: 40px; margin-top: -40px;"),
                             dataTableOutput("overrep"),
                             br(),
                             uiOutput("adapterCon", style="padding-top: 40px; margin-top: -40px;"),
                             imageOutput("adapterContent", height = "600px"),
                             br(),
                             fluidRow(align="center", htmlOutput("topFqc"))
                           )
                         )
                  ))
        ), 
        tabItem(
          tabName = "mappingRates",
          fluidRow(
            uiOutput("sampleMap"),
            column(12,
                   fluidRow(
                     column(6, checkboxGroupInput(
                       "checkGroupQuality", h3("STAR quality metrics"),
                       choices = list("Number of uniquely mapped reads",
                                      "Number of reads mapped to multiple loci",
                                      "Number of reads unmapped"),
                       selected = "Number of uniquely mapped reads")),
                     column(6, align="right", actionButton("mapButton", "Continue"), style= "margin-top:20px;")
                     ),
                   plotOutput("stackedMapping", dblclick = "stackMap_dbl", brush = brushOpts(id = "stackMap_brush", resetOnNew = TRUE), click = "stackMap_click"),
                   br(),
                   fluidRow(column(6, verbatimTextOutput("MapInfoClick")),
                            column(6, verbatimTextOutput("MapInfoBrush"))),
                   br(),
                   fluidRow(align="center", uiOutput("mapSlider")),
                   br(),
                   fluidRow(column(5, align="right", br(), br(), actionButton("calcMap", "Calculated Cutoff")),
                            column(5, align="center", helpText("Download current version of the plot"),
                                   downloadButton("downloadMapping", "Download Plot"))),
                   br(),
                   h3("Uniquely mapped reads"),
                   plotOutput("mappingPlot"),
                   fluidRow(column(12, align= "right", downloadButton("downloadMapBar", "Download plot"))),
                   br(),
                   h3("Metrics of STAR analysis"),
                   dataTableOutput("mappingTab"),
                   fluidRow(align="center", htmlOutput("topMap"))
            ))
        ), tabItem(
          tabName = "geneCounts",
          fluidRow(
            uiOutput("sampleGene"),
            column(12,
                   fluidRow(
                     column(6, checkboxGroupInput(
                       "checkGroupAssigned", h3("FeatureCounts quality metrics"),
                       choices = list("Assigned",
                                      "Unassigned no features",
                                      "Unassigned multimapping",
                                      "Unassigned ambiguity",
                                      "Unassigned unmapped"),
                       selected = "Assigned")),
                     column(6, align="right", actionButton(inputId = "geneButton", label = "Continue"), style= "margin-top:20px;")
                   ),
                   plotOutput("stackedGenePlot", dblclick = "stackGene_dbl", brush = brushOpts(id = "stackGene_brush", resetOnNew = TRUE), click = "stackGene_click"),
                   br(),
                   fluidRow(column(6, verbatimTextOutput("GeneInfoClick")),
                            column(6, verbatimTextOutput("GeneInfoBrush"))),
                   br(),
                   fluidRow(align="center", column(12, uiOutput("geneSlider"))),
                   br(),
                   
                   fluidRow(column(5, align="right", br(), br(), actionButton("calcGene", "Calculated Cutoff")),
                            column(5, align="center", helpText("Download current version of the plot"),
                                   downloadButton("downloadGene", "Download Plot"))),
                   br(),
                   h3("Assigned genes of cell barcodes"),
                   plotOutput("geneCountsPlot"),
                   fluidRow(column(12, align= "right", downloadButton("downloadGeneBar", "Download plot"))),
                   br(),
                   h3("Metrics of featureCounts analysis"),
                   dataTableOutput("geneTab"),
                   fluidRow(align="center", htmlOutput("topGene"))
            ))
        ), tabItem(tabName = "validCells",
                   htmlOutput("validTop"),
                   fluidRow(
                     uiOutput("sampleValid"),
                     column(12,
                            h3("Summary of quality metrics"),
                            plotOutput("umiKnee", dblclick = "kneePlot_dbl", brush = brushOpts(id = "umiKnee_brush", resetOnNew = TRUE), click = "umiKnee_click"),
                            br(),
                            fluidRow(align="center",
                                     column(4, uiOutput("kneeSlider")),
                                     column(2, br(), br(), br(), actionButton("calcCut", "Recommended Cutoff")),
                                     column(6, p("After checking the different metrics of the samples and determining  cells, you can download the expression matrix of the valid cells, for further downstream analysis."),
                                            splitLayout(cellWidths = c("30%", "30%"),
                                              downloadButton("valid", label = "Valid cells"), downloadButton("validParam", label = "Parameters")), style="border:solid; border-color:red; border-radius:4px; padding-bottom:10px; padding-top:10px;")),
                            br(),
                            fluidRow(column(6, verbatimTextOutput("kneeInfoClick")),
                                     column(6, verbatimTextOutput("kneeInfoBrush"))),
                            br(),
                            fluidRow(align="center", helpText("Download current version of the plot")),
                            fluidRow(align="center", downloadButton("downloadKnee", "Download Plot")),
                            br(),
                            h3("Gene counts and mapping metrics of selected cell barcodes"),
                            fluidRow(column(6, plotlyOutput("mappingPie")),
                                     column(6, plotlyOutput("genePie"))),
                            br(),
                            h3("Metrics of UMI counts of all cell barcodes"),
                            dataTableOutput("umiInfoTab"),
                            br(),
                            fluidRow(align="center", htmlOutput("topValid"))
                     ))),
        tabItem(tabName = "contact",
                fluidRow(
                  column(12,
                         strong("Contact:"),
                         p("wasp@computational.bio"),
                         p("andreas.hoek@computational.bio.uni-giessen.de"),
                         br(),
						 strong("Web:"),
                         p(a(href="http://www.computational.bio/software/wasp/", "http://www.computational.bio/software/wasp/")),
						 p(a(href="https://www.uni-giessen.de/fbz/fb08/Inst/bioinformatik", "https://www.uni-giessen.de/fbz/fb08/Inst/bioinformatik")))
                ))
      )
    )
  )


server <- function(input, output, session) {
  
  ### Choosing the uploaded files
  
  # Check if user uploaded new files & call function to verify if all necessary files are uploaded
  observeEvent(input$fileIn,{
    checkInputFiles()
  })
  
  # Verify uploaded files
  checkInputFiles <- reactive({
    if(!is.null(input$fileIn)){
      inputCollapse <- ""
      for(i in input$fileIn$name){
        # Collapse file names for verification
        inputCollapse <- paste(inputCollapse, i, sep = "|")
      }
      if(grepl("_Cell_Numbers.json", inputCollapse) && grepl("_Gene_Counts_Sample.json", inputCollapse) && grepl("_Gene_Counts.json", inputCollapse) && grepl("_Mapping_Rates_Sample.json", inputCollapse) && grepl("_Mapping_Rates.json", inputCollapse) && grepl("_UMI_Counts.json", inputCollapse) && grepl("_fastqc.html", inputCollapse) && grepl("_fastqc.zip", inputCollapse) && grepl("_Demultiplexed.zip", inputCollapse)){
        checkInputFiles <- TRUE
      } else {
        checkInputFiles <- FALSE
      }
    } else {
      checkInputFiles <- FALSE
    }
    checkInputFiles
  })
  
  # Extract name of samples from "_Cell_Numbers.json" file that is included in every valid upload
  getSampleNames <- reactive({
    getSampleNames <- c()
    for (i in seq_along(input$fileIn$name)) {
      if(grepl("_Cell_Numbers.json", input$fileIn$name[i])){
        getSampleNames <- c(getSampleNames, (unique(strsplit(input$fileIn$name[i], "_Cell_Numbers.json"))))
      }
    }
    getSampleNames
  })
  
  # Print sample names for user to website
  output$path <- renderText({
    if(is.null(input$fileIn)){
      paste("No input")}
    else{
      paste(getSampleNames())
    }
  })
  
  
  # Extract sample names
  info <- reactive({
    paste(getSampleNames())
  })
  
  # Print info for user if all necessary files from Snakemake workflow have been uploaded or ask for re-upload
  output$info <- renderPrint({
    if(checkInputFiles() == TRUE){
      print("Correct input files found, continue analysis with menu on the left side.")
    } else if(is.null(input$fileIn) && checkInputFiles() == FALSE){
      print("Please upload input files")
    } else if(!is.null(input$fileIn) && checkInputFiles() == FALSE){
      print("Incorrect input files uploaded, please upload a valid WASP Pre-Processing directory")
    }
  })
  
  
  ## Sidebar Menu
  output$menu <- renderMenu({
    
    if(!checkInputFiles()){
      sidebarMenu(id = "tabs",
                  menuItem("About", tabName = "about"),
                  menuItem("Contact", tabName = "contact")
      )} 
    else if(checkInputFiles()) {
      sidebarMenu(id = "tabs",
                  menuItem("About", tabName = "about"),
                  menuItem("Sample Summaries", tabName = "sampleSummaries"),
                  menuItem("FASTQC Report", tabName = "fastqcResults"),
                  menuItem("Mapping Rates", tabName = "mappingRates"),
                  menuItem("Gene Counts", tabName = "geneCounts"),
                  menuItem("Valid Cells", tabName = "validCells"),
                  menuItem("Contact", tabName = "contact")
      )}
  })
  
  
  # Show information button only on specfic tabs
  observe({
    shinyjs::show("help")
    
    if (is.null(input$tabs) || input$tabs == "contact"){
      shinyjs::hide("help")
    }
  })
  
  
  # If path is defined and all tabs are shown, button to change the active tab should appear
  output$welcomeButton <- renderUI({
    if(checkInputFiles()){
      actionButton(inputId = "welcomeButton", label = "Continue")
    }
  })
  
    
  #If button is pressed, change to Sample Summaries Tab
  observeEvent(input$welcomeButton, {
    updateTabsetPanel(session, "tabs", selected = "sampleSummaries")
  })
  
  
  # Show information for specific tab
  observeEvent(input$help, {
    
    if (is.null(input$tabs)){
      return()
    } 
    
    if (input$tabs == "about"){
      showModal(modalDialog(title = "About",
                            p("For the proper visualization of the data it is important that the preprocessing was accomplished with the snakemake pipeline (sc_pipeline)."),
                            hr(),
                            p("Please select the path to your results directory, that contains the data of the preprocessing. 
                            There should be no other directories in the Results directory except those, the pipeline generated, as these can lead to errors while visualising the data."),
                            br(),
                            p("Example:"),
                            p("../Data/Results/"),
                            easyClose = TRUE)) 
    } else if (input$tabs == "sampleSummaries"){
      showModal(modalDialog(title = "Sample Summaries",
                            p("This tab contains overall information about the samples."),
                            br(),
                            strong("Summary of all samples"),
                            p("Raw reads - The  number of all counted reads in the sample."),
                            p("Reads with valid barcode - The number of reads, that have been validated and were sorted into cell barcodes."),
                            p("Valid percent - The percent of reads, that have been validated."),
                            p("Cell barcodes - The number of barcodes, that have been identified."),
                            hr(),
                            strong("Counted reads of samples"),
                            br(),
                            p("This stacked plot shows the proportion of reads that contain a valid barcode to raw reads in the samples."),
                            hr(),
                            strong("Mapping metrics of STAR analysis for samples"),
                            p("Number of input reads - The number of raw reads in that sample."),
                            p("Average input read length - The average of all read lengths in that sample."),
                            p("Average mapped length - The average mapped read length of all reads in that sample."),
                            p("Mismatch rate per base % - The percentage of possible mismatches per base."),
                            p("Number of uniquely mapped reads - The number of reads that were uniquely mapped to the genome."),
                            p("Uniquely Mapped Reads % - The percentage of reads that were uniquely mapped."),
                            p("Number of reads mapped to multiple loci - The number of reads that were mapped to multiple loci of the genome."),
                            p("% of reads mapped to multiple loci - The percentage of reads that were mapped to multiple loci."),
                            p("Number of reads unmapped - The number of reads that could not be mapped to the genome."),
                            p("% Of Reads Unmapped - The percentage of reads that could not be mapped."),
                            hr(),
                            strong("Mapping metrics of samples"),
                            br(),
                            p("This stacked plot shows the proportion of reads in the samples after mapping."),
                            hr(),
                            strong("Gene count metrics of featureCount analysis for samples"),
                            br(),
                            p("Assigned - The number of reads that were assigned to a specific gene feature."),
                            p("Unassigned ambiguity - The number of reads that could not be assigned to a gene feature due to ambiguity."),
                            p("Unassigned multimapping - The number of reads that were assigned to multiple different gene features."),
                            p("Unassigned no features - The number of reads that were mapped to a region that is not annotated in the annotation file."),
                            p("Unassigned unmapped - The number of reads that could not be mapped to a gene feature."),
                            hr(),
                            strong("Gene count metrics"),
                            br(),
                            p("This stacked plot shows the proportion of reads after counting reads to genomic features."),
                            p("Caution: Due to multimappings read count of featureCounts results can be higher than raw read count."),
                            easyClose = TRUE))
    } else if (input$tabs == "fastqcResults"){
      showModal(modalDialog(title = "FastQC Report",
                            p("The results of the FastQC quality control analysis are shown here."),
                            hr(),
                            p("For further information on the results please visit:"),
                            br(),
                            a(href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/", "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/")
                            ))
      
    } else if (input$tabs == "mappingRates"){
      showModal(modalDialog(title = "STAR Metric Information", 
                            p("This tab shows information gained with STAR mapping."),
                            hr(),
                            strong("Interactive barplot"),
                            p("With checkboxes it is possible to display different metrics in the interactive stacked barplot. Every bar stands for one cell barcode."),
                            p("Through brushing, a specific part can be selected and with double click into the selected part, it is possible to zoom into the plot. To zoom out just double click in the plot (without brushing)."),
                            p("By clicking and brushing in the plot, more information about the cell barcodes can be gained. These will be shown in the text fields under the plot."),
                            strong("Click"),
                            p("Displays the cell barcode and the number of reads at that position."),
                            strong("Brush"),
                            p("Displays the cell barcode range and the number of reads of the selected region."),
                            hr(),
                            strong("Slider"),
                            p("The number of cell barcodes being displayed in the interactive barplot and in the \"uniquely mapped reads\" plot are being controlled with the slider.
                              The initial number of reads being displayed is the calculated number of valid cell barcodes."),
                            hr(),
                            strong("Uniquely mapped reads plot"),
                            p("In this plot the cell barcodes are sorted by the number of their uniquely mapped reads. The mean shows the calculated number of uniquely mapped reads per cell barcode."),
                            hr(),
                            strong("Metrics of STAR analysis for all cell barcodes table"),
                            br(),
                            p("Number of uniquely mapped reads - The number of reads in the cell barcode that were uniquely mapped."),
                            p("Number of reads mapped to multiple loci - The number of reads in that celbarcode that were mapped to multiple loci."),
                            p("Number of reads unmapped - The number of reads in the cell barcode that could not be mapped."),
                            easyClose = TRUE
      ))
    } else if (input$tabs == "geneCounts"){
      showModal(modalDialog(title = "FeatureCounts Metric Information",
                            p("This tab shows information gained with featureCounts analysis."),
                            strong("Interactive barplot"),
                            p("Through selecting the checkboxes it is possible to display different metrics in the interactive stacked barplot. Every bar stands for one cell barcode."),
                            p("Through brushing, a specific part can be selected and with double click it is possible to zoom into this part. To zoom out a double click in the plot is enough 
                              (without brushing)."),
                            p("By clicking and brushing in the plot, more information about the cell barcodes can be gained. These will be shown in the text fields under the plot."),
                            strong("Click"),
                            p("Displays the cell barcode and the number of reads at that position."),
                            strong("Brush"),
                            p("Displays the cell barcode range and the number of reads of the selected region."),
                            hr(),
                            strong("Slider"),
                            p("The number of cell barcodes being displayed in the interactive barplot and in the \"Assigned genes\" plot are being controlled with the slider.
                              The initial number of reads being displayed is the calculated number of valid cell barcodes."),
                            hr(),
                            strong("Assigned genes plot"),
                            p("In this plot the cell barcodes are sorted by the number of their reads assigned to genes. The mean shows the calculated number of assigned reads per cell barcode."),
                            hr(),
                            strong("Metrics of featureCounts analysis for all cell barcodes table"),
                            br(),
                            p("Assigned - The number of reads that were assigned to a specific gene feature in the annotation file."),
                            p("Unassigned No Features - The number of reads that were mapped to a region that is not annotated in the annotation file."),
                            p("Unassigned Multimapping - The number of reads that were assigned to multiple different gene features."),
                            p("Unassigned Ambiguity - The number of reads that could not be assigned to a gene feature due to ambiguity."),
                            easyClose = TRUE))
    } else if (input$tabs == "validCells"){
      showModal(modalDialog(title = "UMI-tools Metric Information",
                            p("This tabs contains information of different analysis steps."),
                            hr(),
                            strong("UMI kneeplot"),
                            p("Here the cell barcodes are displayed as points, they are sorted by the UMI counts in descending order. 
                              A calculated cutoff (red) for valid cell barcodes is displayed. The user is able to determine the cutoff with the slider (blue)."),
                            p("By clicking and brushing in the plot, more information about the cell barcodes can be gained. These will be shown in the text fields under the plot."),
                            strong("Click"),
                            p("Displays the cell barcode and the number of reads at that position."),
                            strong("Brush"),
                            p("Displays the cell barcode range and the number of reads of the selected region."),
                            hr(),
                            strong("Slider"),
                            p("The user has to determine the valid cell barcodes with the slider (blue). We recommend the calculated number of valid cells (red), that can be selected through the \"recommended cutoff\"
                               button. After selecting the number of valid cells with the slider the gene expression matrix can be downloaded as CSV file through the \"Valid cells\" button."),
                            hr(),
                            strong("Pie charts"),
                            p("The mapping rates and gene count metrics of the selected cell barcodes are shown in pie charts."),
                            p("As the pie charts are visualized with the R ploty package, they can be downloaded through a hover element in the plot on the upper right hand side."),
                            hr(),
                            strong("Metrics of UMI counts of all cell barcodes table"),
                            br(),
                            p("UMI Counts - The number of deduplicated UMIs."),
                            p("Genic UMIs - The number of genes in that cell barcode, counted through the deduplicated UMIs."),
                            easyClose = TRUE))
    }
  })
  
  
  ### Sample Summaries Tab
  
  # Anchor to jump to top of tab
  output$sampTop <- renderText({
    "<span id='top'></span>"
  })
  
  
  # Generating lists of JSON files with information of samples
  cell_Numb <- reactive({
    inFiles <- input$fileIn
    cell_Numb <- character()
    for (i in seq_along(inFiles$datapath)) {
      if(grepl("Cell_Numbers.json", inFiles$name[i])){
        cell_Numb <- c(cell_Numb, inFiles$datapath[i])
      }
    }
    cell_Numb
  })
  
  cellMap <- reactive({
    inFiles <- input$fileIn
    cellMap <- character()
    for (i in seq_along(inFiles$datapath)) {
      if(grepl("Mapping_Rates_Sample.json", inFiles$name[i])){
        cellMap <- c(cellMap, inFiles$datapath[i])
      }
    }
    cellMap
  })
  
  cellGene <- reactive({
    inFiles <- input$fileIn
    cellGene <- character()
    for (i in seq_along(inFiles$datapath)) {
      if(grepl("Gene_Counts_Sample.json", inFiles$name[i])){
        cellGene <- c(cellGene, inFiles$datapath[i])
      }
    }
    cellGene
  })
  
  # Sample Summaries Datatable
  output$summaryTab <- DT::renderDataTable({
    
    withProgress(message = "Generating summaries datatable", value = 0, {
      n <- 1
      
      for (j in 1/n){
        
        sampList <- lapply(cell_Numb(), function(i){
          samp <- fromJSON(file = i, method = "C", unexpected.escape = "error", simplify = TRUE)
          rbind(
            c(Sample=samp[["Sample"]],Raw_reads=samp[["Raw_reads"]],Reads_with_valid_barcode=samp[["Reads_with_valid_barcode"]],Valid_percent=percent(as.numeric(samp[["Reads_with_valid_barcode"]]) / as.numeric(samp[["Raw_reads"]])), Cell_barcodes=samp[["Cell_barcodes"]]),
            incProgress(1/n)
          )
        })
      }
      
      sum_Tab <- Reduce(function(x,y) merge(x,y, all=TRUE), sampList)
      colnames(sum_Tab) <- gsub("_", " ", colnames(sum_Tab))
    })
    
    formatCurrency(datatable(sum_Tab, rownames = FALSE, options = list(columnDefs = list(list(className = "dt-right", targets = c(1:4)))))
                   , columns = c(2,3,5), currency = "", interval = 3, mark = ",", digits = 0)

  })
  
  
  # Reactive function to generate bar chart with read counts
  sampSumBar <- reactive({
    
    withProgress(message = "Generating summaries plot", value = 0, {
      n <- 1
      
      for (j in 1/n){
        
        sampList <- lapply(cell_Numb(), function(i){
          samp <- fromJSON(file = i, method = "C", unexpected.escape = "error", simplify = TRUE)
          rbind(
            c(SampleName=samp[["Sample"]],Metrics="Reads without valid barcode",Value=as.numeric(samp[["Raw_reads"]]) - as.numeric(samp[["Reads_with_valid_barcode"]])),
            c(SampleName=samp[["Sample"]],Metrics="Reads with valid barcode",Value=as.numeric(samp[["Reads_with_valid_barcode"]])),
            incProgress(1/n)
          )
        })
      }
      
      
      sum_Tab <- Reduce(function(x,y) rbind(x,y), sampList)
      
      stdf <- as.data.frame(sum_Tab,stringsAsFactors = c(T,T,F))
      stdf$Value <- as.numeric(as.character(stdf$Value))
    })
    
    ggplot(stdf, aes(x= SampleName, y=Value, fill = Metrics)) + geom_bar(stat = "identity", position = "fill") + ggtitle("Counted reads of samples") + scale_fill_manual(values = c("dodgerblue2", "gold2")) + theme_classic() +
      labs(y="Percent") + theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), axis.title = element_text(size = 13), legend.title=element_text(size=13), legend.text=element_text(size=12), legend.position = "top") + scale_y_continuous(labels = scales::percent_format())+
      geom_text(aes(label = format(Value, big.mark = ",")), position = position_fill(vjust = 0.5))
  })
  
  
  # Bar chart for with read counts
  output$sumBarChart <- renderPlot({
    sampSumBar()
  })
  
  
  # Download bar chart with read counts
  output$downloadSampSum <- downloadHandler(
    
    filename = function(){
      paste("CountedReadsBarChart", "pdf", sep=".")
    },
    content = function(file){
      withProgress(message = "Processing", value = 0, {
        n <- 1
        
        for (j in 1/n){
          pdf(file, width = 15, height = 10)
          plot(sampSumBar())
          dev.off()
          incProgress(1/n)
        }
      })
    }
  )
  
  # Generating datatable with STAR results
  output$mappingSum <- DT::renderDataTable({
    
    withProgress(message = "Generating mapping datatable", value = 0, {
      n <- 1
      
      for (j in 1/n){
        
        sampList <- lapply(cellMap(), function(i){
          samp <- fromJSON(file = i, method = "C", unexpected.escape = "error", simplify = TRUE)
          samp <- as.data.frame(samp)
        })
        incProgress(1/n)
      } 
      
      mapSum <- Reduce(function(x,y) merge(x,y, all = TRUE), sampList)
      datatable(mapSum)
      
      colnames(mapSum) <- gsub("(X.|[.])", "%", colnames(mapSum))
      colnames(mapSum) <- gsub("_", " ", colnames(mapSum))
      
    })
    
    formatCurrency(datatable(mapSum[,c("Sample", "Number of input reads", "Average input read length", "Average mapped length", "Mismatch rate per base %", "Number of uniquely mapped reads", "Uniquely mapped reads %", "Number of reads mapped to multiple loci", "% of reads mapped to multiple loci", "Number of reads unmapped", "% of reads unmapped")],
                             rownames = FALSE, options = list(
                               columnDefs = list(list(className = "dt-right", targets = c(2:10)))
                             )), columns = c(2, 3, 6, 8, 10), currency = "", interval = 3, mark = ",", digits = 0)
  })
  
  
  # Generating bar chart with STAR results of all samples
  sampMapBar <- reactive({
    
    withProgress(message = "Generating bar chart (mapping)", value = 0, {
      n <- 1
      
      for (j in 1/n){

        sampList <- lapply(cellMap(), function(i){
          samp <- fromJSON(file = i, method = "C", unexpected.escape = "error", simplify = TRUE)
          rbind(
            c(SampleName=samp[["Sample"]],Metrics="Number of uniquely mapped reads",Value=as.numeric(samp[["Number_of_uniquely_mapped_reads"]])),
            c(SampleName=samp[["Sample"]],Metrics="Number of reads mapped to multiple loci",Value=as.numeric(samp[["Number_of_reads_mapped_to_multiple_loci"]])),
            c(SampleName=samp[["Sample"]],Metrics="Number of reads unmapped", Value=as.numeric(samp[["Number_of_reads_unmapped"]])),
            incProgress(1/n)
          )
        })
      }
      
      sum_Tab <- Reduce(function(x,y) rbind(x,y), sampList)
      
      stdf <- as.data.frame(sum_Tab,stringsAsFactors = c(T,T,F))
      stdf$Value <- as.numeric(as.character(stdf$Value))
      
      ggplot(stdf, aes(x= SampleName, y=Value, fill = Metrics)) + geom_bar(stat = "identity", position = "fill") + ggtitle("Mapping metrics of samples") + scale_fill_manual(values = c("dodgerblue2", "gold2", "seagreen3")) + theme_classic() +
        labs(y="Percent") + theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), axis.title = element_text(size = 13), legend.title=element_text(size=13), legend.text=element_text(size=12), legend.position = "top") + scale_y_continuous(labels = scales::percent_format())+
        geom_text(aes(label = format(Value, big.mark = ",")), position = position_fill(vjust = 0.5))
    }) 
  })
  
  
  # Stacked bar chart with STAR results
  output$mapSumChart <- renderPlot({
    sampMapBar()
  })
  
  
  # Download bar chart of STAR results for all samples
  output$downloadSampMap <- downloadHandler(
    filename = function(){
      paste("MappingSummaryBarChart", "pdf", sep=".")
    },
    content = function(file){
      withProgress(message = "Processing", value = 0, {
        n <- 1
        
        for (j in 1/n){
          pdf(file, width = 15, height = 10)
          plot(sampMapBar())
          dev.off()
          incProgress(1/n)
        }
      })
    }
  )
  
  
  # Generating datatable with featureCounts results
  output$geneSum <- DT::renderDataTable({
    
    withProgress(message = "Generating gene datatable", value = 0, {
      n <- 1
      
      for (j in 1/n){

        sampList <- lapply(cellGene(), function(i){
          samp <- fromJSON(file = i, method = "C", unexpected.escape = "error", simplify = TRUE)
          samp <- as.data.frame(samp)
        })
        incProgress(1/n)
      }
      
      geneSum <- Reduce(function(x,y) merge(x,y, all = TRUE), sampList)
      names(geneSum) <- gsub("_", " ", names(geneSum))
    })
    
    formatCurrency(datatable(geneSum[,c("Sample", "Assigned", "Unassigned ambiguity", "Unassigned multimapping", "Unassigned no features", "Unassigned unmapped")], 
                             rownames = FALSE), columns = c(2:6), currency = "", interval = 3, mark = ",", digits = 0)
  })  
  
  
  # Generating stacked bar chart for featureCounts metrics
  sampGeneBar <- reactive({
    
    withProgress(message = "Generating bar chart (gene)", value = 0, {
      n <- 1
      
      for (j in 1/n){
        
        sampList <- lapply(cellGene(), function(i){
          samp <- fromJSON(file = i, method = "C", unexpected.escape = "error", simplify = TRUE)
          rbind(
            c(SampleName=samp[["Sample"]],Metrics="Assigned",Value=as.numeric(samp[["Assigned"]])),
            c(SampleName=samp[["Sample"]],Metrics="Unassigned ambiguity",Value=as.numeric(samp[["Unassigned_ambiguity"]])),
            c(SampleName=samp[["Sample"]],Metrics="Unassigned multimapping", Value=as.numeric(samp[["Unassigned_multimapping"]])),
            c(SampleName=samp[["Sample"]],Metrics="Unassigned no features", Value=as.numeric(samp[["Unassigned_no_features"]])),
            c(SampleName=samp[["Sample"]],Metrics="Unassigned unmapped", Value=as.numeric(samp[["Unassigned_unmapped"]])),
            incProgress(1/n)
          )
        })
      }
      
      sum_Tab <- Reduce(function(x,y) rbind(x,y), sampList)
      
      stdf <- as.data.frame(sum_Tab,stringsAsFactors = c(T,T,F))
      stdf$Value <- as.numeric(as.character(stdf$Value))
      
      ggplot(stdf, aes(x= SampleName, y=Value, fill = Metrics)) + geom_bar(stat = "identity", position = "fill") + ggtitle("Gene count metrics of samples") + scale_fill_manual(values = c("dodgerblue2", "gold2", "seagreen3", "mediumorchid2", "#bc3c4d")) + theme_classic() +
        labs(y="Percent") + theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), axis.title = element_text(size = 13), legend.title=element_text(size=13), legend.text=element_text(size=12), legend.position = "top") + scale_y_continuous(labels = scales::percent_format())+
        geom_text(aes(label = format(Value, big.mark = ",")), position = position_fill(vjust = 0.5))
    })
  }) 
  
  
  # Stacked bar chart for featureCounts metrics
  output$geneSumChart <- renderPlot({
    sampGeneBar()
  })
  
  
  # Download bar chart of featureCounts results for all samples
  output$downloadSampGene <- downloadHandler(
    filename = function(){
      paste("GeneSummaryBarChart", "pdf", sep=".")
    },
    content = function(file){
      withProgress(message = "Processing", value = 0, {
        n <- 1
        
        for (j in 1/n){
          
          pdf(file, width = 15, height = 10)
          plot(sampGeneBar())
          dev.off()
          incProgress(1/n)
        }
      })
    }
  )
  
  
  # Adding a go to top button
  output$topSum <- renderText({
    "<hr><span id='anchorSum'></span>
    <a href='#top' class = 'btn btn-default' role='button'>Go to top</a>"
  })
  
  
  # Button to change from sample summaries tab to fastqc results tab
  observeEvent(input$sumButton, {
    updateTabsetPanel(session, "tabs", selected = "fastqcResults")
  })
  
  
  
  ## FASTQC Results
  
  # Generating tabnames for every sample
  output$fcSamples <- renderUI({
    do.call(tabsetPanel, c(id = "fc", lapply(info(), function(i){
      tabPanel(title = paste0(i))
    })))
  })
  
  
  # Print tabname that is active
  observeEvent(input$fc, {
    print(input$fc)
    
    # Unzip FastQC files
    inFiles <- input$fileIn
    
    withProgress(message = "Loading data", value = 0, {
      n <- 1
      for (j in 1/n){
        for (i in seq_along(inFiles$datapath)) {
          if(grepl(input$fc, inFiles$name[i]) && grepl("_fastqc.zip", inFiles$name[i])){
            unzip(zipfile = inFiles$datapath[i], exdir = paste0("fastqc_results/", input$fc, "/FastQC_Files"))
          }
          incProgress(1/n)
        }
      }
    })
  })  
  
  
  # FastQC Sidebar
  
  # FastQC sidebar: Basic Statistics
  # Link to jump to basic statistics with icon that reflects the quality
  output$basic <- renderUI({
    if (!is.null(input$fc) || length(input$fc) > 0) {
      a(href="#basicSta", fastqcIcon(fastqcSum(input$fc)[1]), "Basic Statistics")
    }
  })
  
  
  # FastQC sidebar: Per base sequence quality
  # Link to jump to per base sequence quality with icon that reflects the quality
  output$baseSeq <- renderUI({
    if (!is.null(input$fc) || length(input$fc) > 0) {
      a(href="#perBaseQua", fastqcIcon(fastqcSum(input$fc)[2]),"Per base sequence quality")
    }
  })
  
  
  # FastQC sidebar: Per tile sequence quality
  # Link to jump to per tile sequence quality with icon that reflects the quality
  output$tileQua <- renderUI({
    if (!is.null(input$fc) || length(input$fc) > 0) {
      a(href="#tileQuali", fastqcIcon(fastqcSum(input$fc)[3]),"Per tile sequence quality")
    }
  })
  
  
  # FastQC sidebar: Per sequence quality scores
  # Link to jump to per sequence quality scores with icon that reflects the quality
  output$seqQua <- renderUI({
    if (!is.null(input$fc) || length(input$fc) > 0) {
      a(href="#qualitySco", fastqcIcon(fastqcSum(input$fc)[4]),"Per sequence quality scores")
    }
  })
  
  
  # FastQC sidebar: Per base sequence content
  # Link to jump to per base sequence content with icon that reflects the quality
  output$seqCon <- renderUI({
    if (!is.null(input$fc) || length(input$fc) > 0) {
      a(href="#sequenceCont", fastqcIcon(fastqcSum(input$fc)[5]),"Per base sequence content")
    }
  })
  
  
  # FastQC sidebar: Per sequence GC content
  # Link to jump to per sequemce GC content with icon that reflects the quality
  output$gcCon <- renderUI({
    if (!is.null(input$fc) || length(input$fc) > 0) {
      a(href="#gcCont", fastqcIcon(fastqcSum(input$fc)[6]),"Per sequence GC content")
    }
  })
  
  
  # FastQC sidebar: Per base N content
  # Link to jump to per base N content with icon that reflects the quality
  output$nCon <- renderUI({
    if (!is.null(input$fc) || length(input$fc) > 0) {
      a(href= "#nCont", fastqcIcon(fastqcSum(input$fc)[7]),"Per base N content")
    }
  })
  
  
  # FastQC sidebar: Sequence Length Distribution
  # Link to jump to sequence length distribution with icon that reflects the quality
  output$lenDist <- renderUI({
    if (!is.null(input$fc) || length(input$fc) > 0) {
      a(href ="#lengthDist", fastqcIcon(fastqcSum(input$fc)[8]),"Sequence Length Distribution")
    }
  })
  
  
  # FastQC sidebar: Sequence Duplication Levels
  # Link to jump to sequence duplication levels with icon that reflects the quality
  output$dupLev <- renderUI({
    if (!is.null(input$fc) || length(input$fc) > 0) {
      a(href="#duplicationLev", fastqcIcon(fastqcSum(input$fc)[9]),"Sequence Duplication Levels")
    }
  })
  
  
  # FastQC sidebar: Overrepresented sequences
  # Link to jump to overrepresented sequences with icon that reflects the quality
  output$overrepSeq <- renderUI({
    if (!is.null(input$fc) || length(input$fc) > 0) {
      a(href="#overrepresentedSeq", fastqcIcon(fastqcSum(input$fc)[10]),"Overrepresented sequences")
    }
  })
  
  
  # FastQC sidebar: Adapter Content
  # Link to jump to adapter content with icon that reflects the quality
  output$adapCon <- renderUI({
    if (!is.null(input$fc) || length(input$fc) > 0) {
      a(href="#adapterCon", fastqcIcon(fastqcSum(input$fc)[11]),"Adapter Content")
    }
  })
  
  
  # FastQC Main Panel
  
  # FastQC main: Basic Statistics
  # Anchor for basic statistics
  output$basicSta <- renderUI({ h3("Basic Statistics") })
  
  # Generating datatable with FastQC results of basic statistics
  output$basicstats <- DT::renderDataTable({
    
    withProgress(message = "Generating datatable", value = 0, {
      n <- 1
      
      for (j in 1/n){
        statsTab <- data.frame()
        
        if(!is.null(input$fc) & length(input$fc) > 0){
          fastqc_lines <- read_lines(paste0("fastqc_results/", input$fc, "/FastQC_Files/", input$fc, "_fastqc/fastqc_data.txt", sep = ""))[4:10]
          
          data <- lapply(fastqc_lines, function(i){
            l <- unlist(strsplit(i, split = "\t"))
            rbind(c(l))
          })
          
          statsTab <- Reduce(function(x,y) rbind(x,y), data)
          colnames(statsTab) <- c("Measure", "Value")
          statsTab[4:7,2] <- format(as.numeric(statsTab[4:7,2]), big.mark = ",")
          incProgress(1/n)
        }
      }
    }) 
    
    datatable(statsTab)
    
  }) 
  
  
  # FastQC main: Per Base Sequence Quality Image
  # Anchor for per base sequence quality image
  output$perBaseQua <- renderUI({h3("Per base sequence quality")})
  
  # FastQC results for per base sequence quality image
  output$perBaseQuality <- renderImage({
    return(list(src = (paste0("fastqc_results/", input$fc, "/FastQC_Files/", input$fc, "_fastqc/Images/per_base_quality.png")), contentType = "Image/png", alt="Quality per base"))
  }, deleteFile = FALSE)
  
  
  # FastQC main: Per tile sequence quality
  # Anchor for per tile sequence quality
  output$tileQuali <- renderUI({h3("Per tile sequence quality")})
  
  # FastQC results for per tile sequence quality
  output$tileQuality <- renderImage({
    return(list(src = (paste0("fastqc_results/", input$fc, "/FastQC_Files/", input$fc, "_fastqc/Images/per_tile_quality.png")), contentType = "Image/png", alt="Quality per tile"))
  }, deleteFile = FALSE)
  
  
  # FastQC main: Per sequence quality scores
  # Anchor for per sequence quality scores
  output$qualitySco <- renderUI({h3("Per sequence quality scores")})
  
  # FastQC results for per sequence quality scores
  output$qualityScores <- renderImage({
    return(list(src = (paste0("fastqc_results/", input$fc, "/FastQC_Files/", input$fc, "_fastqc/Images/per_sequence_quality.png")), contentType = "Image/png", alt="Quality per sequence"))
  }, deleteFile = FALSE)
  
  
  # FastQC main: Per base sequence content
  # Anchor for per base sequence content
  output$sequenceCont <- renderUI({h3("Per base sequence content")})
  
  # FastQC results for per base sequence content
  output$sequenceContent <- renderImage({
    return(list(src = (paste0("fastqc_results/", input$fc, "/FastQC_Files/", input$fc, "_fastqc/Images/per_base_sequence_content.png")), alt="Per base sequence content"))
  }, deleteFile = FALSE)
  
  
  # FastQC main: Per sequence GC content
  # Anchor for per sequence GC content
  output$gcCont <- renderUI({h3("Per base GC content")})
  
  # FastQC Results for per sequence GC content
  output$gcContent <- renderImage({
    return(list(src = (paste0("fastqc_results/", input$fc, "/FastQC_Files/", input$fc, "_fastqc/Images/per_sequence_gc_content.png")), alt="Per sequence GC content"))
  }, deleteFile = FALSE)
  
  
  # FastQC main: Per base N content
  # Anchor for per base n content
  output$nCont <- renderUI({h3("Per base N content")})
  
  # FastQC results for per base n content
  output$nContent <- renderImage({
    return(list(src = (paste0("fastqc_results/", input$fc, "/FastQC_Files/", input$fc, "_fastqc/Images/per_base_n_content.png")), alt="Per base N content"))
  }, deleteFile = FALSE)
  
  
  # FastQC main: Sequence Length Distribution
  # Anchor for sequence length distribution
  output$lengthDist <- renderUI({h3("Sequence Length Distribution")})
  
  # FastQC results for sequence length distribution
  output$lengthDistribution <- renderImage({
    return(list(src = (paste0("fastqc_results/", input$fc, "/FastQC_Files/", input$fc, "_fastqc/Images/sequence_length_distribution.png")), alt="Sequence Length Distribution"))
  }, deleteFile = FALSE)
  
  
  # FastQC main: Sequence Duplication Levels
  # Anchor to sequence duplication levels
  output$duplicationLev <- renderUI({h3("Sequence Duplication Levels")})
  
  # FastQC results for sequence duplication levels
  output$duplication <- renderImage({
    return(list(src = (paste0("fastqc_results/", input$fc, "/FastQC_Files/", input$fc, "_fastqc/Images/duplication_levels.png")), alt="Sequence Duplication Levels"))
  }, deleteFile = FALSE)
  
  
  # FastQC main: Overrepresented Sequences
  # Anchor to overrepresented sequences
  output$overrepresentedSeq <- renderUI({h3("Overrepresented sequences")})
  
  # Generating datatable with FastQC results for overrepresented sequences
  output$overrep <- DT::renderDataTable({
    
    if(!is.null(input$fc) || length(input$fc) > 0){
      
      withProgress(message = "Generating table", value = 0, {
        n <- 1
        
        for (j in 1/n){
          fastqc_lines <- read_lines(paste0("fastqc_results/", input$fc, "/FastQC_Files/", input$fc, "_fastqc/fastqc_data.txt", sep = ""))
          metrics <- grep(">>", fastqc_lines)
          overrepLines <- fastqc_lines[(metrics[19]+2):(metrics[20]-1)]
          
          data <- lapply(overrepLines, function(i){
            l <- unlist(strsplit(i, split = "\t"))
            rbind(c(l))
          })
          
          incProgress(1/n)
          
          statsTab <- Reduce(function(x,y) rbind(x,y), data)
          colnames(statsTab) <- c("Sequence", "Count", "Percentage", "Possible Source")
        }
      })
      
      formatCurrency(datatable(statsTab, options = list(
        columnDefs = list(list(className = "dt-right", targets = 2))),
        rownames = FALSE), columns = 2, currency = "", interval = 3, mark = ",", digits = 0)
    }
  }) 
  
  
  # FastQC main: Adapter Content
  # Anchor to adapter content
  output$adapterCon <- renderUI({h3("Adapter Content")})
  
  # FastQC results for adapter content
  output$adapterContent <- renderImage({
    return(list(src = (paste0("fastqc_results/", input$fc, "/FastQC_Files/", input$fc, "_fastqc/Images/adapter_content.png")), alt="Adapter Content"))
  }, deleteFile = FALSE)
  
  
  # Button to change from fastqc results tab to mapping rates tab
  observeEvent(input$fastButton, {
    updateTabsetPanel(session, "tabs", selected = "mappingRates")
  })
  
  
  # Adding a go to top button
  output$topFqc <- renderText({
    "<hr><span id='anchordFqc'></span>
    <a href='#fcSamples' class = 'btn btn-default' role='button'>Go to top</a>"
  })
  
  
  
  # Reactive values for stacked bar charts
  stackedPlots <- reactiveValues()
  
  
  
  # STAR Mapping Results
  
  # Generating tabnames for every sample
  output$sampleMap <- renderUI({
    do.call(tabsetPanel, c(id = "map", lapply(info(), function(i){
      tabPanel(title = paste0(i))
    })))
  })
  
  
  # Print tabname that is active
  observeEvent(input$map, {
    print(input$map)
  })
  
  
  # Reading in mapping information of sample
  sampleMapData <- reactive({
    inFiles <- input$fileIn
    mappingFile <- NULL 
    for (i in seq_along(inFiles$datapath)) {
      if(grepl(input$map, inFiles$name[i]) && grepl("_Mapping_Rates.json", inFiles$name[i])){
        mappingFile <- inFiles$datapath[i]
      }
    }
    fromJSON(file = mappingFile, method = "C", unexpected.escape = "error", simplify = TRUE)
  })
  
  
  # Reactive function for reading-in data and generating a matrix for the stacked bar chart 
  sortedMap <- reactive(if(!is.null(input$map) || length(input$map) > 0){
    
    withProgress(message = "Reading data", value = 0, {
      n <- 1
      
      for (j in 1/n){
        sampleTab <- as.data.frame(sampleMapData(), row.names("Cell_barcode"))
        sampleTabSorted <- sampleTab[order(-sampleTab$Number_of_uniquely_mapped_reads),]
        names(sampleTabSorted)<- gsub("_", " ", names(sampleTabSorted))
        incProgress(1/n)
      }
      })
    
    return(sampleTabSorted)
  }) 
  
  
  # Reactive values for brushing event in stacked bar charts
  rangesMap <- reactiveValues(x = NULL, y = NULL)
  
  
  # Ranges for brushing event in stacked bar chart
  observeEvent(input$stackMap_dbl, {
    
    brush <- input$stackMap_brush
    
    if (!is.null(brush)) {
      rangesMap$x <- c(brush$xmin, brush$xmax)
      rangesMap$y <- c(brush$ymin, brush$ymax)
    } else {
      rangesMap$x <- NULL
      rangesMap$y <- NULL
    }
    
  })

  
  # Generating stacked bar chart for STAR metrics
  output$stackedMapping <- renderPlot({
    
    if(is.null(sortedMap())){return()}
    if(!length(input$mapSlider) > 0){return()}
    if(!is.null(input$map) || length(input$map) > 0){
      
      withProgress(message = "Generating stacked plot", value = 0, {
        n <- 1
        
        for (j in 1/n){
          samp2 <- sortedMap()[,-1]
          rownames(samp2) <- sortedMap()[,1]
          
          stackedList <- lapply(rownames(samp2)[1:input$mapSlider], function(i){
            Reduce(function(x,y) rbind(x,y),
                   lapply(input$checkGroupQuality, function(j){
                     rbind(c(as.character(i), j, samp2[i,j]))
                   }))
          })
          incProgress(1/n)
          
          stackedTab <- Reduce(function(x,y) rbind(x,y), stackedList)
        }
        
        if(!is.null(input$checkGroupQuality)){
          colnames(stackedTab) <- c("Cellbarcode", "Metrics", "Value")
          
          stdf <- as.data.frame(stackedTab,stringsAsFactors = c(T,T,F))
          stdf$Value <- as.numeric(as.character(stdf$Value))
          
          gm <- ggplot(stdf, aes(x= reorder(Cellbarcode, -Value), y=Value, fill=Metrics, width=.75)) + geom_bar(stat = "identity") +
            ggtitle("Mapping quality metrics") + coord_cartesian(xlim = rangesMap$x, ylim = rangesMap$y, expand = FALSE) + scale_fill_manual(values = c("#bc3c4d", "#3c8dbc", "#bc3c8d")) + labs(x= "Cell barcodes", y="Number of reads") + theme_classic() +
            theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), axis.text.x = element_blank(), legend.title=element_text(size=13), legend.text=element_text(size=12), legend.position = "top", axis.title = element_text(size = 13))
          
          stackedPlots$gm <- gm
          gm
        }
      })
    }
  })
  
  
  # Info box for click event in stacked bar chart
  output$MapInfoClick <- renderText({
    xy_cell <- function(e) {
      if(is.null(e)) return("\n")
      paste0("Cell barcode range: ", (sortedMap()[,1])[round(e$x, 0)], "\nNumber of reads: ", format(round(e$y, 2), big.mark = ","), "\n")
    }
    paste0(
      "Click - Current Selection:\n", xy_cell(input$stackMap_click)
    )
  })
  
  
  # Info box for brush event in stacked bar chart
  output$MapInfoBrush <- renderText({
    
    range_xy_cell <- function(e) {
      if(is.null(e)) return("\n")
      
      paste0("Cell barcode range:\nFrom: ", (sortedMap()[,1])[round(e$xmin, 0)], "\nTo: ", (sortedMap()[,1])[round(e$xmax, 0)], "\t\n",
             "Number of reads - min: ", format(round(e$ymin, 2), big.mark = ","), "\nNumber of reads - max: ", format(round(e$ymax, 2), big.mark = ","))
    }
    
    paste0("Brush - Current Selection:\n", range_xy_cell(input$stackMap_brush))
  })
  
  
  # Mapping Rates slider to control the number of shown cell barcodes in the stacked bar chart
  # Initial value is one and a half times the calculated number of valid cells
  output$mapSlider <- renderUI({
    
    if(length(input$map) > 0){
      
      umiTab <- kneePlot(input$map)
      
      rawdiff <- diff(umiTab$log_lib_size)/diff(umiTab$barcode_rank)
      inflection <- which(rawdiff == min(rawdiff[100:length(rawdiff)], na.rm = TRUE))
      valid <- round(inflection + (inflection*0.5), digits = 0)
      
      sliderInput("mapSlider", h3("Number of barcodes", br(), h5("Warning: Selecting a high value might lead to an increased calculation time.")), min=1, max = nrow(sortedMap()), value = valid, step = 1, width = "700px")
      
    }
  })
  
  
  # Updating slider in stacked bar chart, if button for recommended cut off is used
  observeEvent(input$calcMap, {
    
    umiTab <- kneePlot(input$map)
    
    withProgress(message = "Updating slider", value = 0, {
      n <- 1
      
      for (j in 1/n){
        rawdiff <- diff(umiTab$log_lib_size)/diff(umiTab$barcode_rank)
        inflection <- which(rawdiff == min(rawdiff[100:length(rawdiff)], na.rm = TRUE))
        incProgress(1/n)
      }
    })
    
    updateSliderInput(session, "mapSlider", value = inflection)
  })
  
  
  # Download stacked bar chart
  output$downloadMapping <- downloadHandler(
    
    filename = function(){
      paste(input$map, "StackedMapping.pdf", sep="_")
    },
    content = function(aFile){
      
      withProgress(message = "Processing", value = 0, {
        n <- 1
        
        for (j in 1/n){
          pdf(aFile, width = 15, height = 10)
          print(stackedPlots$gm)
          dev.off()
          incProgress(1/n)
        }
      })
    }
  )
  
  
  # Reactive function to generate histogram
  uniquelyMap <- reactive({
    
    if(!length(input$mapSlider) > 0){return()}
    
    sampleTab <- sortedMap()
    
    ggplot(sampleTab[1:input$mapSlider,], aes(sampleTab[,2][1:input$mapSlider])) + geom_histogram(bins=length(sampleTab[,2][1:input$mapSlider]), fill="dodgerblue2", col="black") + theme_classic() + 
      geom_vline(aes(xintercept = mean(sampleTab[,2][1:input$mapSlider]), color = paste(format(round(mean(sampleTab[,2][1:input$mapSlider]), digits = 2), big.mark = ","), " Reads per cell barcode"))) + scale_color_manual(name = "Mean", values = "red") + 
      ggtitle(paste("Uniquely mapped reads of cell barcodes of", input$map)) + labs(x="Uniquely mapped reads", y="Frequency") + 
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), legend.title=element_text(size=13), legend.text=element_text(size=12), axis.title = element_text(size = 13))
  })
  
  
  # Histogram for uniquely mapped reads of selected cell barcodes
  output$mappingPlot <- renderPlot({
    
    if(!is.null(input$map) || length(input$map) > 0){
      uniquelyMap()
    }
  }) 
  
  
  # Download for histogram of uniquely mapped reads
  output$downloadMapBar <- downloadHandler(
    filename = function(){
      paste("UniquelyMapped", "pdf", sep=".")
    },
    content = function(file){
      
      withProgress(message = "Processing", value = 0, {
        n <- 1
        
        for (j in 1/n){
          pdf(file, width = 15, height = 10)
          plot(uniquelyMap())
          dev.off()
          incProgress(1/n)
        }
      })
    })
  
  
  # Generating datatable with STAR results
  output$mappingTab <- DT::renderDataTable({
    
    if(!is.null(input$map) || length(input$map) > 0){
      
      withProgress(message = "Generating datatable", value = 0, {
        n <- 1
        
        for (j in 1/n){
          
          sampleTab <- as.data.frame(sampleMapData(), row.names("Cell_barcode"))
          names(sampleTab) <- gsub("_", " ", names(sampleTab))
          incProgress(1/n)
        }
        
        formatCurrency(datatable(sampleTab[, c("Cell barcode", "Number of uniquely mapped reads", "Number of reads mapped to multiple loci", "Number of reads unmapped")], 
                                 rownames = FALSE, options = list(
                                   order = list(list(1, "desc")))), columns = c(2:4), currency = "", interval = 3, mark = ",", digits = 0)
      })
    } 
  })
  
  
  # Adding a go to top button
  output$topMap <- renderText({
    "<hr><span id='anchorMap'></span>
    <a href='#sampleMap' class = 'btn btn-default' role='button'>Go to top</a>"
  })
  
  
  observeEvent(input$mapButton, {
    updateTabsetPanel(session, "tabs", selected = "geneCounts")
  })
  
  
  
  # Results from featureCounts Analysis
  
  # Generating tabnames for every sample
  output$sampleGene <- renderUI({
    do.call(tabsetPanel, 
            c(id = "gene", lapply(info(), function(i){
              tabPanel(title = paste0(i))
            })
            ))
  })
  
  
  # Print tabname that is active
  observeEvent(input$gene, {
    print(input$gene)
  })
  
  
  # Reading in gene counts information of sample
  sampleGeneData <- reactive({
    inFiles <- input$fileIn
    geneFile <- NULL 
    for (i in seq_along(inFiles$datapath)) {
      if(grepl(input$gene, inFiles$name[i]) && grepl("_Gene_Counts.json", inFiles$name[i])){
        geneFile <- inFiles$datapath[i]
      }
    }
    fromJSON(file = geneFile, method = "C", unexpected.escape = "error", simplify = TRUE)
  })

  
  # Reactive function for reading-in data and generating a matrix for the stacked bar chart
  sortedGene <- reactive(if(!is.null(input$gene) || length(input$gene) > 0){
    
    withProgress(message = "Generating stacked plot", value = 0, {
      n <- 1
      
      for (j in 1/n){
        sampleTab <- as.data.frame(sampleGeneData(), row.names("Cellbarcode"))
        sampleTabSorted <- sampleTab[order(-sampleTab$Assigned),]
        names(sampleTabSorted)<- gsub("_", " ", names(sampleTabSorted))
        incProgress(1/n)
      }
    })
    
    return(sampleTabSorted)
  })
  
  
  # Reactive values for the brushing event in the stacked bar chart
  rangesGene <- reactiveValues(x = NULL, y = NULL)
  
  
  # Ranges for brushing event in stacked bar chart
  observeEvent(input$stackGene_dbl, {
    
    brush <- input$stackGene_brush
    
    if (!is.null(brush)) {
      rangesGene$x <- c(brush$xmin, brush$xmax)
      rangesGene$y <- c(brush$ymin, brush$ymax)
    } else {
      rangesGene$x <- NULL
      rangesGene$y <- NULL
    }
    
  })
  
  
  # Generating stacked bar chart for featureCounts metrics
  output$stackedGenePlot <- renderPlot({
    
    if(is.null(sortedGene())){return()}
    
    if(!length(input$geneSlider) > 0){return()}
    
    if(!is.null(input$gene) || length(input$gene) > 0){
      
      withProgress(message = "Generating stacked plot", value = 0, {
        n <- 1
        
        for (j in 1/n){
          samp2 <- sortedGene()[,-1]
          rownames(samp2) <- sortedGene()[,1]
          
          stackedList <- lapply(rownames(samp2)[1:input$geneSlider], function(i){
            
            Reduce(function(x,y) rbind(x,y),
                   
                   lapply(input$checkGroupAssigned, function(j){
                     rbind(c(as.character(i), j, samp2[i,j]))
                   })
            )
          })
          incProgress(1/n)
          
          stackedTab <- Reduce(function(x,y) rbind(x,y), stackedList)
        }
        
        if(!is.null(input$checkGroupAssigned)){
          colnames(stackedTab) <- c("Cellbarcode", "Metrics", "Value")
          
          stdf <- as.data.frame(stackedTab,stringsAsFactors = c(T,T,F))
          stdf$Value <- as.numeric(as.character(stdf$Value))
          
          gg <- ggplot(stdf, aes(x= reorder(Cellbarcode, -Value), y=Value, fill = Metrics, width=0.75)) + ggtitle("Counted gene quality metrics")  + geom_bar(stat = "identity") + 
            coord_cartesian(xlim = rangesGene$x, ylim = rangesGene$y, expand = FALSE) + scale_fill_manual(values = c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc")) + labs(x= "Cell barcodes", y="Number of reads") + theme_classic() +
            theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), axis.text.x = element_blank(), legend.title=element_text(size=13), legend.text=element_text(size=12), legend.position = "top", axis.title = element_text(size = 13))
          
          stackedPlots$gg <- gg
          gg
        }
      })
    }
  })
  
  
  # Info box for click event in stacked bar chart
  output$GeneInfoClick <- renderText({
    
    xy_cell <- function(e) {
      if(is.null(e)) return("\n")
      
      paste0("Cell barcode range: ", (sortedGene()[,1])[round(e$x, 0)], "\nNumber of reads: ", format(round(e$y, 2), big.mark = ","), "\n")
    }
    
    paste0("Click - Current Selection:\n", xy_cell(input$stackGene_click))
  })
  
  
  # Info box for brush event in stacked bar chart
  output$GeneInfoBrush <- renderText({
    
    range_xy_cell <- function(e) {
      if(is.null(e)) return("\n")
      
      paste0("Cell barcode range:\nFrom: ", (sortedGene()[,1])[round(e$xmin, 0)], "\nTo: ", (sortedGene()[,1])[round(e$xmax, 0)], "\t\n",
             "Number of reads - min: ", format(round(e$ymin, 2), big.mark = ","), "\nNumber of reads - max: ", format(round(e$ymax, 2), big.mark = ","))
    }
    
    paste0("Brush - Current Selection:\n", range_xy_cell(input$stackGene_brush))
  })
  
  
  # Gene Counts slider to control the number of shown cell barcodes in the stacked bar chart
  # Initial value is one and a half times the calculated number of valid cells
  output$geneSlider <- renderUI({
    
    if(length(input$gene) > 0){
      
      umiTab <- kneePlot(input$gene)
      
      rawdiff <- diff(umiTab$log_lib_size)/diff(umiTab$barcode_rank)
      inflection <- which(rawdiff == min(rawdiff[100:length(rawdiff)], na.rm = TRUE))
      valid <- round(inflection + (inflection*0.5), digits = 0)
      
      sliderInput("geneSlider", h3("Number of barcodes", br(), h5("Warning: Selecting a high value might lead to an increased calculation time.")), min=1, max = nrow(sortedGene()), value = valid, step = 1, width = "700px")
    }
  })
  
  
  # Updating slider in stacked bar chart, if button for recommended cut off is used
  observeEvent(input$calcGene, {
    
    umiTab <- kneePlot(input$gene)
    
    withProgress(message = "Updating slider", value = 0, {
      n <- 1
      
      for (j in 1/n){
        rawdiff <- diff(umiTab$log_lib_size)/diff(umiTab$barcode_rank)
        inflection <- which(rawdiff == min(rawdiff[100:length(rawdiff)], na.rm = TRUE))
        incProgress(1/n)
      }
    })
    updateSliderInput(session, "geneSlider", value = inflection)
  })
  
  
  # Download stacked bar chart of featureCounts metrics
  output$downloadGene <- downloadHandler(
    
    filename = function(){
      paste(input$gene, "StackedGene.pdf", sep="_")
    },
    content = function(aFile){
      
      withProgress(message = "Processing", value = 0, {
        n <- 1
        
        for (j in 1/n){
          pdf(aFile, width = 15, height = 10)
          print(stackedPlots$gg)
          dev.off()
          incProgress(1/n)
        }
        })
    }
  )
  
  
  # Reactive function for generating histogram of assigned reads
  assignedGene <- reactive({
    
    if(!length(input$geneSlider) > 0){return()}
    
    sampleTab <- sortedGene()
    
    ggplot(sampleTab[1:input$geneSlider,], aes(sampleTab[,2][1:input$geneSlider])) + geom_histogram(bins=length(sampleTab[,2][1:input$geneSlider]), fill="dodgerblue2", col="black") + theme_classic() +
      geom_vline(aes(xintercept = mean(sampleTab[,2][1:input$geneSlider]), color = paste(format(round(mean(sampleTab[,2][1:input$geneSlider]), digits = 2), big.mark = ","), " Reads per cell barcode"))) + scale_color_manual(name = "Mean", values = "red") +
      ggtitle(paste("Assigned reads of cell barcodes of", input$gene)) + labs(x="Assigned reads", y="Frequency") +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), legend.title=element_text(size=13), legend.text=element_text(size=12), axis.title = element_text(size = 13))
  })
  
  
  # Histogram for assigned reads in selected cell barcodes 
  output$geneCountsPlot <- renderPlot({
    
    if(!is.null(input$gene) || length(input$gene) > 0){
      assignedGene() 
    }
  })
  
  
  # Download histogram of assigned reads of cell barcodes
  output$downloadGeneBar <- downloadHandler(
    filename = function(){
      paste("UniquelyAssigned", "pdf", sep=".")
    },
    content = function(file){
      
      withProgress(message = "Processing", value = 0, {
        n <- 1
        
        for (j in 1/n){
          pdf(file, width = 15, height = 10)
          plot(assignedGene())
          dev.off()
          incProgress(1/n)
        }
      })
    })
  
  
  # Generating datatable with featureCounts results
  output$geneTab <- DT::renderDataTable({
    
    if(!is.null(input$gene) || length(input$gene) > 0){
      
      withProgress(message = "Generating datatable", value = 0, {
        n <- 1
        
        for (j in 1/n){
          sampleTab <- as.data.frame(sampleGeneData(), row.names("Cellbarcodes"))
          
          names(sampleTab) <- gsub("_", " ", names(sampleTab))
          incProgress(1/n)
        }
        
        formatCurrency(datatable(sampleTab[, c("Cell barcode", "Assigned", "Unassigned multimapping", "Unassigned ambiguity", "Unassigned no features", "Unassigned unmapped")],
                                 rownames = FALSE, options = list(
                                   order = list(list(1, "desc")))), columns = c(2:6), currency = "", interval = 3, mark = ",", digits = 0)
      })
    }
  })
  
  
  # Button to change from Gene Counts tab to Valid Vells tab
  observeEvent(input$geneButton, {
    updateTabsetPanel(session, "tabs", selected = "validCells")
  })
  
  
  # Adding a go to top button
  output$topGene <- renderText({
    "<hr><span id='anchorGene'></span>
    <a href='#sampleGene' class = 'btn btn-default' role='button'>Go to top</a>"
  })
  
  
  
  ## Valid Cells Results
  
  # Generating tabnames for every sample
  output$sampleValid <- renderUI({
    do.call(tabsetPanel, c(id = "umi", lapply(info(), function(i){
      tabPanel(title = paste0(i))
    })))
  })
  
  
  # Print tabname that is active
  observeEvent(input$umi, {
    print(input$umi)
  })
  
  
  # Anchor to jump to top of tab
  output$validTop <- renderText({
    "<span id='vtop'></span>"
  })
  
  
  # Generating UMI knee plot
  rangesKnee <- reactiveValues(x = NULL, y = NULL)
  
  output$umiKnee <- renderPlot({
    
    if(!is.null(input$umi) || length(input$umi) > 0){
      umiKneeFunc()
    }
  })
  
  
  # Generating knee slider for knee plot
  output$kneeSlider <- renderUI({
    
    sliderInput("umiKneeSlider", h3("Valid cells cutoff", br(), h5("Warning: Selecting a high value might lead to an increased calculation time.")), min=1, max = nrow(kneePlot(input$umi)), value = 400, step = 1)
  })
  
  
  # Updating slider in UMI knee plot, if button for recommended cut off is used
  observeEvent(input$calcCut, {
    
    umiTab <- kneePlot(input$umi)
    
    withProgress(message = "Updating slider", value = 0, {
      n <- 1
      
      for (j in 1/n){
        rawdiff <- diff(umiTab$log_lib_size)/diff(umiTab$barcode_rank)
        inflection <- which(rawdiff == min(rawdiff[100:length(rawdiff)], na.rm = TRUE))
        incProgress(1/n)
      }
    })
    
    updateSliderInput(session, "umiKneeSlider", value = inflection)
  })
  
  
  # Download button for gene expressionmatrix
  df <- data.frame(Gene=character())
  
  output$valid <- downloadHandler(
    
    
    
    filename = function(){
      paste(input$umi, "GeneMatrix.csv", sep = "_")
    },
    content = function(aFile){
      
      withProgress(message = "Processing", value = 0, {
        n <- 1
        
        inFiles <- input$fileIn
        # Unzip demultiplex files (TSV file for each detected barcode)
        for (i in seq_along(inFiles$datapath)) {
          if(grepl(input$umi, inFiles$name[i]) && grepl("_Demultiplexed.zip", inFiles$name[i])){
            unzip(zipfile = inFiles$datapath[i], exdir = paste0("demultiplexed/", input$umi))
          }
        }
        for (j in 1/n){
          filenames <- row.names(kneePlot(input$umi))[1:input$umiKneeSlider]
          
          splitFileNames <- strsplit(filenames, "_Zelle_")
          samp <- unlist(splitFileNames)[2*(1:length(filenames))-1]
          cellBc <- unlist(splitFileNames)[2*(1:length(filenames))]
          
          for (i in cellBc){
            test.df <- read.table(file = paste("demultiplexed/", input$umi, "/Temporary/", input$umi, "/Demultiplexed_", input$umi, "/", substr(i, 1, 2), "/", substr(i, 3, 4), "/", input$umi, "_Zelle_", i, ".tsv", sep = ""), header = TRUE, sep = "\t")
            incProgress(1/n)
            Gene <- test.df[,1]
            Value <- test.df[,2]
            dft <- data.frame(Gene, Value)
            dft$Gene <- as.character(dft$Gene)
            colnames(dft)[2] <- paste0(input$umi, "_", i)
            df <- merge(df, dft, by="Gene", all = TRUE)
          }

          df[is.na(df)] <- 0
        }
        
        write.table(df, aFile, row.names = FALSE, sep = ",")
      })
    })
  
  
  # Download textfile with parameters for valid cell barcodes
  output$validParam <- downloadHandler(
    
    filename = function(){
      paste(input$umi, "Parmeters.txt", sep = "_")
    },
    content = function(aFile){
      
      withProgress(message = "Processing", value = 0, {
        n <- 1
        
        snakefileFound = FALSE
        
        inFiles <- input$fileIn
        # Unzip demultiplex files (TSV file for each detected barcode)
        for (i in seq_along(inFiles$datapath)) {
          if(grepl("Snakefile", inFiles$name[i])){
            snakefile <- read_lines(inFiles$datapath[i])
            snakefileFound = TRUE
          }
        }
        
        for (j in 1/n){
          samp <- paste("Sample:", input$umi)
          #gen <- paste("Genome:", list.files(path = paste(dirname(path()), "/Reference/", sep = ""), pattern = ".fa"))
          if(snakefileFound){
            fqc <- paste("FastQC:", trimws(snakefile[63]))
            star <- paste("STAR:", trimws(snakefile[91]))
            fc <- paste("featureCounts:", trimws(snakefile[108]))
            ut <- paste("UMI-tools:", trimws(snakefile[202]))
          } else {
            fqc <- paste("FastQC:", "Please Upload Snakefile for these parameters")
            star <- paste("STAR:", "Please Upload Snakefile for these parameters")
            fc <- paste("featureCounts:", "Please Upload Snakefile for these parameters")
            ut <- paste("UMI-tools:", "Please Upload Snakefile for these parameters")
          }
          vbc <- paste("Number of selected barcodes: ", input$umiKneeSlider)
          incProgress(1/n)
        }
        
        writeLines(c(samp, fqc, star, fc, ut, vbc), aFile) #gen,
      })
    })
  
  
  # Info box for click event in UMI knee plot
  output$kneeInfoClick <- renderText({
    
    xy_str <- function(e) {
      if(is.null(e)) return("\n")
      paste0("Cell barcode range: ", (row.names(kneePlot(input$umi))[round(e$x, 0)]), "\n(Log) Genic UMI count: ", round(e$y, 1), "\n")
    }
    paste0(
      "Click - Current Selection: \n", xy_str(input$umiKnee_click)
    )
  })
  
  
  # Ranges for brushing event in UMI knee plot
  observeEvent(input$kneePlot_dbl, {
    
    brush <- input$umiKnee_brush
    
    if (!is.null(brush)) {
      rangesKnee$x <- c(brush$xmin, brush$xmax)
      rangesKnee$y <- c(brush$ymin, brush$ymax)
    } else {
      rangesKnee$x <- NULL
      rangesKnee$y <- NULL
    }
    
  })
  
  
  # Info box for brush event in UMI knee plot
  output$kneeInfoBrush <- renderText({
    
    range_xy_cell <- function(e) {
      if(is.null(e)) return("\n")
      
      paste0("Cell barcode range:\nFrom: " ,(row.names(kneePlot(input$umi))[round(e$xmin, 0)]), "\nTo: ", (row.names(kneePlot(input$umi))[round(e$xmax, 0)]), "\t\n",
             "(Log) Genic UMI count - min: ", round(e$ymin, 2), "\n(Log) Genic UMI count - max: ", round(e$ymax, 2))
    }
    paste0("Brush - Current Selection:\n", range_xy_cell(input$umiKnee_brush))
  })
  
  
  # Download UMI knee plot
  output$downloadKnee <- downloadHandler(
    
    filename = function(){
      paste(input$umi, "ValidCells.pdf", sep="_")
    },
    content = function(aFile){
      withProgress(message = "Processing", value = 0, {
        n <- 1
        
        for (j in 1/n){
          ggsave(aFile, plot = umiKneeFunc(), device = "pdf", width = 15, height = 10)
        }
        incProgress(1/n)
      })
    }
  )
  
  
  # Generating pie chart for mapping information of selected cell barcodes
  output$mappingPie <- renderPlotly({
    
    if(!is.null(input$umi) || length(input$umi) > 0){
      
      withProgress(message = "Generating pie plot (mapping)", value = 0, {
        n <- 1
        
        for(j in 1/n){
          
          inFiles <- input$fileIn
          mappingRatesFile <- NULL 
          for (i in seq_along(inFiles$datapath)) {
            if(grepl(input$umi, inFiles$name[i]) && grepl("_Mapping_Rates.json", inFiles$name[i])){
              mappingRatesFile <- inFiles$datapath[i]
            }
          }
          
          sampleTab <- as.data.frame(fromJSON(file = mappingRatesFile, method = "C", unexpected.escape = "error", simplify = TRUE))
          
          row.names(sampleTab) <- sampleTab$Cell_barcode
          
          umiTab <- kneePlot(input$umi)
          
          stackedList <- lapply(row.names(umiTab)[1:input$umiKneeSlider], function(i){
            rbind(
              c(Cellbarcode=as.character(i), Metric="Uniquely mapped reads", Value=as.numeric(sampleTab[i,"Number_of_uniquely_mapped_reads"])),
              c(Cellbarcode=as.character(i), Metric="Reads mapped to multiple loci", Value=as.numeric(sampleTab[i,"Number_of_reads_mapped_to_multiple_loci"])),
              c(Cellbarcode=as.character(i), Metric="Unmapped reads", Value=as.numeric(sampleTab[i,"Number_of_reads_unmapped"]))
            )
          })
          
          stackedTab <- Reduce(function(x,y) rbind(x,y), stackedList)
          
          incProgress(1/n)
          
          stdf <- as.data.frame(stackedTab,stringsAsFactors = c(T,T,F))
          stdf$Value <- as.numeric(as.character(stdf$Value))
        }
        
        plot_ly(stdf, labels = stdf$Metric, values = stdf$Value, type = "pie") %>%
          
          layout(title= "Mapping rates in valid cell barcodes")
      })
    }
  })
  
  
  # Generating pie chart for gene counts information of selected cell barcodes
  output$genePie <- renderPlotly({
    
    if(!is.null(input$umi) || length(input$umi) > 0){
      
      withProgress(message = "Generating pie plot (gene counts)", value = 0, {
        n <- 1
        
        for(j in 1/n){      
          
          inFiles <- input$fileIn
          geneCountsFile <- NULL 
          for (i in seq_along(inFiles$datapath)) {
            if(grepl(input$umi, inFiles$name[i]) && grepl("_Gene_Counts.json", inFiles$name[i])){
              geneCountsFile <- inFiles$datapath[i]
            }
          }
          
          sampleTab <- as.data.frame(fromJSON(file = geneCountsFile, method = "C", unexpected.escape = "error", simplify = TRUE))
          
          row.names(sampleTab) <- sampleTab$Cell_barcode
          
          umiTab <- kneePlot(input$umi)
          
          stackedList <- lapply(row.names(umiTab)[input$umiKneeSlider], function(i){
            rbind(
              c(Cellbarcode=as.character(i), Metric="Assigned", Value=as.numeric(sampleTab[i,"Assigned"])),
              c(Cellbarcode=as.character(i), Metric="Unassigned no features", Value=as.numeric(sampleTab[i,"Unassigned_no_features"])),
              c(Cellbarcode=as.character(i), Metric="Unassigned multimapping", Value=as.numeric(sampleTab[i,"Unassigned_multimapping"])),
              c(Cellbarcode=as.character(i), Metric="Unassigned ambiguity", Value=as.numeric(sampleTab[i,"Unassigned_ambiguity"])),
              c(Cellbarcode=as.character(i), Metric="Unassigned unmapped", Value=as.numeric(sampleTab[i,"Unassigned_unmapped"]))
            )
          })
          
          stackedTab <- Reduce(function(x,y) rbind(x,y), stackedList)
          
          incProgress(1/n)
          
          stdf <- as.data.frame(stackedTab,stringsAsFactors = c(T,T,F))
          stdf$Value <- as.numeric(as.character(stdf$Value))
        }
        
        plot_ly(stdf, labels = stdf$Metric, values = stdf$Value, type = "pie") %>%
          layout(title= "Counted genes in valid cell barcodes")
      })
    }
  })
  
  
  # Generating datatable for UMI counts
  output$umiInfoTab <- DT::renderDataTable({
    
    if(!is.null(input$umi) || length(input$umi) > 0){
      
      withProgress(message = "Generating datatable", value = 0, {
        n <- 1
        
        for (j in 1/n){      
          
          inFiles <- input$fileIn
          umiCountsFile <- NULL 
          for (i in seq_along(inFiles$datapath)) {
            if(grepl(input$umi, inFiles$name[i]) && grepl("_UMI_Counts.json", inFiles$name[i])){
              umiCountsFile <- inFiles$datapath[i]
            }
          }
          
          sampleData <- fromJSON(file = umiCountsFile, method = "C", unexpected.escape = "error", simplify = TRUE)
          sampleTab <- as.data.frame(sampleData, row.names("Cell_barcode"))
          
          incProgress(1/n)
          
          names(sampleTab) <- gsub("_", " ", names(sampleTab))
          datatable(sampleTab)
        }
        
        formatCurrency(datatable(sampleTab[, c("Cell barcode", "UMI counts", "Gene counts")], rownames = FALSE,
                                 options =list(order = list(list(1, "decs")))), columns = c(2:3),currency = "", interval = 3, mark = ",", digits = 0)
      })
    }
  })
  
  
  # Adding a go to top button
  output$topValid <- renderText({
    "<hr><span id='anchorValid'></span>
    <a href='#sampleValid' class = 'btn btn-default' role='button'>Go to top</a>"
  })
  

  # Function for calculating the valid cells
  kneePlot <- function(sample){
    
    if(!is.null(sample) || length(sample) > 0){
      
      withProgress(message = "Calculating valid cell cutoff", value = 0, {
        n <- 1
        
        for (j in 1/n){
          
          inFiles <- input$fileIn
          umiCountsFile <- NULL 
          for (i in seq_along(inFiles$datapath)) {
            if(grepl(sample, inFiles$name[i]) && grepl("_UMI_Counts.json", inFiles$name[i])){
              umiCountsFile <- inFiles$datapath[i]
            }
          }
          
          sampleUmi <- fromJSON(file = umiCountsFile, method = "C", unexpected.escape = "error", simplify = TRUE)
          sampleTab <- as.data.frame(sampleUmi)
          
          barcode_rank <- rank(-sampleTab$UMI_counts)
          log_lib_size <- log10(sampleTab$UMI_counts)
          
          o <- order(barcode_rank)
          log_lib_size <- log_lib_size[o]
          barcode_rank <- barcode_rank[o]
          
          incProgress(1/n)
          
          umiTab <- data.frame(barcode_rank, log_lib_size)
          row.names(umiTab) <- sampleTab$Cell_barcode[o]
        }
      })
      return(umiTab)
    }} 
  
  # Function
  
  # Function for reading-in summary of FASTQC report
  fastqcSum <- function(fc){
    
    if (!is.null(input$fc) & length(input$fc) > 0) {
      
      summ <- read.csv(file = paste0("fastqc_results/", input$fc, "/FastQC_Files/", fc, "_fastqc/summary.txt"), sep = "\t", header = FALSE, col.names = c("Quality", "Metric", "Sample"))
      return(summ$Quality)
    }
  }
  
  # Function for icon selection that reflects the quality of the sample statistics
  fastqcIcon <- function(metric){
    
    if (metric == "WARN"){
      icon("fas fa-exclamation-circle")
    }
    else if (metric == "PASS"){
      icon("fas fa-check-circle")
    }
    else {
      icon("fas fa-times-circle")
    }
  }
  
  # Function to generate UMI knee plot
  umiKneeFunc <- function(){
    umiTab <- kneePlot(input$umi)
    
    withProgress(message = "Generating knee plot", value = 0, {
      n <- 1
      
      for (j in 1/n){
        rawdiff <- diff(umiTab$log_lib_size)/diff(umiTab$barcode_rank)
        inflection <- which(rawdiff == min(rawdiff[100:length(rawdiff)], na.rm = TRUE))
        incProgress(1/n)
      }
    })
    
    ggplot(umiTab, aes(x=barcode_rank, y=log_lib_size)) + geom_point() + ggtitle("UMIs per cell barcode") +
      coord_cartesian(xlim = rangesKnee$x, ylim = rangesKnee$y) +
      labs(x= "Cell barcodes in descending genic UMI count order", y="(log) Genic UMI count") + theme_classic() +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), axis.title = element_text(size = 13)) + 
      geom_vline(aes(xintercept = inflection, colour = paste0(inflection))) + scale_color_manual(name = "Valid cells cutoff\n(calculated)", values = "red") + 
      geom_vline(xintercept = input$umiKneeSlider, color = "blue") + scale_x_continuous(expand = c(0, 0.1)) + scale_y_continuous(expand=expand_scale(mult=c(0.1,0.1)))
  }
  
  # Remove unzipped directories containing FastQC results & demultiplexed TSV files
  session$onSessionEnded(function(){unlink(c("demultiplexed", "fastqc_results"), recursive = TRUE, force = TRUE)})
  # Terminate session if browser is closed
  session$onSessionEnded(stopApp)
  
}

shinyApp(ui, server)