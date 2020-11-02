# Shiny
library(shiny)
library(shinydashboard)

# analysis
library(dplyr)
library(DT)
library(ggplot2)
library(plotly)
library(scater)
library(Seurat)
library(SingleCellExperiment)

# size of uploadable data increased to 500MB
options(shiny.maxRequestSize=500*1024^2)

ui <- dashboardPage(
  dashboardHeader(title="WASP"),
  
  dashboardSidebar(
    # HTML / CSS code
    tags$head(
      tags$style(
        HTML(
          "#summary{
            margin-top: 200px;
          }
          #goToUpload{
            display: block;
            margin-left: auto;
            font-size: 16px;
          }
          #goToNorm{
            display: block; 
            margin-left: auto; 
          }
          .glyphicon.glyphicon-forward{
            color: #DC143C;
          }
          .glyphicon.glyphicon-info-sign{
            color: #3c8dbc;
          }
          #downloadAllMarkers{
            display: block;
            margin-left: auto;
            margin-right: auto;
            width: 50%;
          } 
          #downloadSpecMarkers{
            display: block;
            margin-left: auto;
            margin-right: auto;
            width: 50%;
          }
          #showInfoUpload{
            display: block;
            margin-left: auto;
            margin-right: auto;
            float: left;
          }
          .norm_button{
            display: block;
            margin-left:auto;
            float: right;
          }
          .clust_button{
            display: block;
            margin-left:auto;
            float: right;
          }
          .marker_button{
            display: block;
            margin-left:auto;
            float: right;
          }
          #showInfoResults{
            display: block;
            margin-left: auto;
          }
          .modal-dialog{
            width: 1000px;
          }
          .modal-header{
            background-color: #3c8dbc;
            color: #fff;
          }
          table {
            margin: 15px;
          }
          table#custom, table#custom th, table#custom td {
            border-collapse: collapse;
          }
          table#custom th, table#custom td {
            padding: 5px;
            text-align: left;
          }
          table#custom tr:nth-child(even) {
            background-color: #dbe3ec;
          }
          table#custom tr:nth-child(odd) {
            background-color: #ecf0f5;
          }
          table#custom th {
            background-color: #3c8dbc;
            color: white;
          }"
        )
      )
    ),
    sidebarMenu(
      id="sidebar",
      menuItem("Welcome", tabName = "welcome", icon = icon("home", lib="glyphicon")),
      menuItem("Data Upload", tabName = "upload", icon = icon("file", lib="glyphicon")),
      sidebarMenuOutput(outputId="filter_tab"),
      sidebarMenuOutput(outputId="norm_tab"),
      sidebarMenuOutput(outputId="cluster_tab"),
      sidebarMenuOutput(outputId="marker_tab"),
      sidebarMenuOutput(outputId = "result_tab"),
      menuItem("Impressum", tabName="impressum", icon = icon("copyright-mark", lib="glyphicon"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName="welcome",
              h1("Welcome,", style="color: #3c8dbc;"),
              div(style="font-size: 16px;", br(), "and thank you for using", span(style="color: #3c8dbc;", strong("WASP")), "(Web-accessible single-cell platform).", br(), br(),
              "The purpose of this application is the visualization of processed single-cell RNA sequencing data and includes the following steps:",
              p(style="margin-left: 20px;", "- filtering low-quality cells and genes", br(),
              "- normalizing the data", br(),
              "- clustering of cells", br(),
              "- finding possible markergenes for the estimated clusters"), 
              "Every step outputs plots, which are possible to be downloaded on your local computer, to visualize the outcome.", br(), br(),
              "To start with the analysis and upload your data, press the button below."),
              br(),
              actionButton("goToUpload", "Start Analysis", icon = icon("forward", lib="glyphicon"))
      ),
      tabItem(tabName="upload",
              fluidRow(column(4, fileInput("exprmatrix", "Please select your .csv file containing the expression matrix you want to analyze",
                                           buttonLabel = "Browse...", placeholder = "No file selected",
                                           accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
                       column(8, htmlOutput("expr_explanation"))),
              fluidRow(column(4, fileInput("annotable", "Please select your .csv file containing the annotation table if available (optional)",
                                           buttonLabel = "Browse...", placeholder = "No file selected",
                                           accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                              uiOutput('ui.annoDropdown')),
                       column(8, htmlOutput("anno_explanation"))),
              fluidRow(column(12, htmlOutput("expr_upload_ok"))),
              fluidRow(column(12, htmlOutput("anno_upload_ok"))),
              br(),
              fluidRow(column(1, uiOutput('ui.info'))),
              br(),
              fluidRow(column(2, uiOutput('ui.radio'))),
              fluidRow(column(1, uiOutput('ui.action')))      
      ),
      tabItem(tabName="results",
              fluidRow(
                column(10, uiOutput("ui.downloadButtons")),
                column(2, actionButton("showInfoResults", "Get Information", icon=icon("info-sign", lib="glyphicon")))
              ),
              br(),
              fluidRow(
                column(9, uiOutput("ui.tabsetResults")),
                column(3, tableOutput("summary")))
      ),
      tabItem(tabName="filter",
              fluidRow(
                column(6, downloadButton("filter_man", "Download Plots")),
                column(6, uiOutput("ui.startNorm"))),
              br(),
              fluidRow(
                column(6, uiOutput("ui.sliderUmi")),
                column(6, uiOutput("ui.sliderGene"))),
              fluidRow(
                column(6, plotOutput("umihisto", height="400px")),
                column(6, plotOutput("genehisto", height="400px"))),
              br(),
              textOutput("cellsremain"),
              br(),
              h3("Gene Filter"),
              fluidRow(
                column(6, uiOutput("ui.numTranscripts")), 
                column(6, uiOutput("ui.numCells"))),
              fluidRow(column(12, plotOutput("genedistrib", height="500px"))),
              br(),
              fluidRow(column(6, textOutput("genesremain")))
      ),
      tabItem(tabName="normalize",
              fluidRow(
                column(5, downloadButton("norm_man", "Download Plots"),
                       downloadButton("normtable_man", "Download Normalization Table")),
                uiOutput("ui.startCluster")
              ),
              br(),
              fluidRow(
                column(6, uiOutput("ui.significantPCs")),
                column(6, uiOutput("ui.resolution"))
              ),
              tabsetPanel(type="tabs",
                          tabPanel("Elbow Plot", plotOutput("elbow_plot", height="650px")),
                          tabPanel("Jackstraw Plot", plotOutput("jackstraw_plot", height="650px")),
                          tabPanel("VizPCA Plot", plotOutput("vizpca_plot", height="650px")),
                          tabPanel("PCHeatmap", plotOutput("pcheatmap", height="650px")))
      ),
      tabItem(tabName="cluster",
              fluidRow(
                column(5, downloadButton("cluster_man", "Download Plots")),
                uiOutput("ui.startMarker")
              ),
              br(),
              fluidRow(
                column(6, uiOutput("ui.findMarkers")),
                column(6, uiOutput("ui.compareMarkers"))),
              checkboxInput("allMarkers", label = "Find markergenes for all clusters", value = TRUE),
              fluidRow( 
                column(6, numericInput("numPval", label="p-value", value=0.01, min=0.01, max=1.00, step=0.01)),
                column(6, numericInput("numLogfc", label="log2 fold change", value=0.25, min=0.01, max=1.00, step=0.01))),
              h3("Clustered cells"),
              checkboxInput("scaleCells_clust1", label = "Scale cells by number of genes", value = TRUE),
              tabsetPanel(type="tabs",
                          tabPanel("2D",
                                   tabsetPanel(type="tabs",
                                               tabPanel("t-SNE of the clustered cells",
                                                        fluidRow(column(12, plotOutput("tsne_clustered", height="550px")))),
                                               tabPanel("UMAP of the clustered cells",
                                                        fluidRow(column(12, plotOutput("umap_clustered", height="550px")))),
                                               tabPanel("PCA of the clustered cells",
                                                        fluidRow(column(12, plotOutput("pca_clustered", height="550px")))),
                                               tabPanel("MDS of the clustered cells",
                                                        fluidRow(column(12, plotOutput("mds_clustered", height="550px")))))),
                          tabPanel("3D",
                                   tabsetPanel(type="tabs",
                                               tabPanel("t-SNE",
                                                        fluidRow(column(12, plotlyOutput("tsne_clust3d", height="550px")))),
                                               tabPanel("UMAP",
                                                        fluidRow(column(12, plotlyOutput("umap_clust3d", height="550px")))),
                                               tabPanel("PCA",
                                                        fluidRow(column(12, plotlyOutput("pca_clust3d", height="550px")))),
                                               tabPanel("MDS",
                                                        fluidRow(column(12, plotlyOutput("mds_clust3d", height="550px"))))))
                         )
      ),
      tabItem(tabName="marker",
              fluidRow(
                column(3, downloadButton("marker_man", "Download Plots")),
                column(6, uiOutput("ui.downloadMarkersAll"),
                       uiOutput("ui.downloadMarkersSpec")),
                column(3, tags$div(class="marker_button", actionButton("showInfoMark", "Get Information", icon=icon("info-sign", lib="glyphicon")),
                                   actionButton("goToResult", "View All Results", icon = icon("forward", lib="glyphicon"))))
              ),
              br(),
              uiOutput("ui.tabsetMarkers")
      ),
      tabItem(tabName="impressum",
              h2("Impressum"),
              fluidRow(column(12, htmlOutput("imp")))
      )
    )
  )
)