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
      menuItem("Contact", tabName="contact", icon = icon("copyright-mark", lib="glyphicon"))
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
      tabItem(tabName="contact",
              h2("Contact:"),
              fluidRow(column(12, htmlOutput("imp")))
      )
    )
  )
)

server <- function(input, output, session) {
  
  ### preparation for data analysis
  
  ## Switching to upload tab
  observeEvent(input$goToUpload, {
    updateTabItems(session, "sidebar", selected="upload")
  })
  
  ## Explanation of the expression matrix
  output$expr_explanation <- renderText({
    paste0("<p style='padding-left:15px;'>The expression matrix must contain cell IDs (column names) as well as gene IDs (row names).</p>",
           "<table id='custom' style='width:40%;'>
             <tr>
               <th></th>
               <th>cell_ID1</th>
               <th>cell_ID2</th>
               <th>...</th>
             </tr>
             <tr>
               <td><b>gene_ID1</b></td>
               <td>14</td>
               <td>3</td>
               <td>...</td>
             </tr>
             <tr>
               <td><b>gene_ID2</b></td>
               <td>6</td>
               <td>32</td>
               <td>...</td>
             </tr>
             <tr>
               <td><b>...</b></td>
               <td>...</td>
               <td>...</td>
               <td>...</td>
             </tr>
            </table>")
  })
  
  ## Explanation of the annotation table
  output$anno_explanation <- renderText({
    paste0("<table id='custom' style='width:27%; float:left;'>
             <tr>
               <th>cell_IDs</th>
               <th>celltype</th>
             </tr>
             <tr>
               <td>cell_ID1</td>
               <td>metadata</td>
             </tr>
             <tr>
               <td>cell_ID2</td>
               <td>metadata</td>
             </tr>
             <tr>
               <td>...</td>
               <td>...</td>
             </tr>
            </table>",
           "<table style='width:60%; float:left;'>
             <tr>
               <td style='line-height: 1.7'>The annotation table usually contains metadata of the individual cells</td>
             </tr>
             <tr>
               <td style='line-height: 1.7'>and only needs to be uploaded if the information are important for the analysis.</td>
             </tr>
             <tr>
               <td style='line-height: 1.7'>In the following steps the metadata will be used to distinguish the cells presented</td>
             </tr>
             <tr>
               <td style='line-height: 1.7'>in different plots (e.g. t-SNE or PCA plots). Every step of the analysis will be</td>
             </tr>
             <tr>
               <td style='line-height: 1.7'>carried out even without providing an annotation table.</td>
             </tr>
            </table>"
    )
  })

  ## Information Button Upload
  output$ui.info <- renderUI({
    if(is.null(uploadFile())) {
      return()
    }
    actionButton("showInfoUpload", "Get Information", icon=icon("info-sign", lib="glyphicon"))
  })
  
  observeEvent(input$showInfoUpload, {
    showModal(modalDialog(
      title=HTML("<b>Type of Analysis:</b>"),
      HTML("<div style='line-height: 1.7'><b>To start the analysis you have to decide between the two provided types:</b><br>
           - the automatic analysis, where your data will be processed automatically<br> 
             with certain default parameters<br>
           - or the manual analysis, where you will be led through every step of the data processing<br>
             and can adjust the parameters to better fit your samples.<br><br>
           The automatic analysis is a quick but complete run through the workflow, the manual analysis is recommended for experienced users who wish for cutomization.</div>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  ## Information Button Normalization
  observeEvent(input$showInfoNorm, {
    showModal(modalDialog(
      title=HTML("<b>Additional information regarding Normalization:</b>"),
      HTML("<div style='line-height: 1.7'><h4>Significant principal components</h4>
            'To overcome the extensive technical noise in any single gene for scRNA-seq data, [cells are clustered] based on their PCA scores, 
            with each PC essentially representing a ‘metagene’ that combines information across a correlated gene set.'<sup>1</sup><br>
            'We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.'<sup>2</sup><br>
            The default parameter has been calculated with the data shown in the elbow plot. You as a user can decide whether to trust this calculation
            or readjust the significant principal components based on the jackstraw plot, elbow plot and the PC Heatmap displayed below.<br>
            <h4>Resolution for clustering</h4>
            '[The] resolution parameter [...] sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. 
            We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.'<sup>1</sup>
            <h4>Elbow plot</h4>
            'A ranking of principle components based on the percentage of variance explained by each one.'<sup>1</sup>
            The elbow in this plot indicates the last significant principal component.
            <h4>Jackstraw plot</h4>
            'A visualization tool for comparing the distribution of p-values for each [principal component] with a uniform distribution (dashed line).'<sup>1</sup>
            <h4>VizPCA</h4>
            'Visualize top genes associated with reduction components.'<sup>3</sup> (principal components)             
            <h4>PC Heatmap</h4>
            'Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.'<sup>3</sup><br>
            '[It also] allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses.'<sup>1</sup>
            <br><br>
            <sup>1</sup> <a href='https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html' target='_blank'>Seurat - Guided Clustering Tutorial V3.0</a> <br>
            <sup>2</sup> <a href='https://satijalab.org/seurat/v1.4/pbmc3k_tutorial.html' target='_blank'>Seurat - Guided Clustering Tutorial V1.4</a> <br>
            <sup>3</sup> <a href='https://cran.r-project.org/web/packages/Seurat/Seurat.pdf' target='_blank'>Seurat Package Manual</a>
           </div>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  ## Information Button Clustering
  observeEvent(input$showInfoClust, {
    showModal(modalDialog(
      title=HTML("<b>Additional information:</b>"),
      HTML("<div style='line-height: 1.7'><h4>Cluster for marker gene analysis</h4>
            A specific cluster can be selected in the dropdown menu for which the respective marker genes will be calculated. For this calculation the selection of the second dropdown menu will be taken into account. 
            The gene expression of both clusters will be looked at and compared to find significant marker genes specifically for the one cluster that has been selected.<br><br>
            There is also the additional possibility to calculate the marker genes of every cluster of the analysis. 
            <h4>p-value</h4>
            Simply put, in this case the p-value describes the probability of a gene not being a significant marker gene. The lower the p-value the better. 
            The threshold you can set will define which markers will be looked at in the analysis: 'Only return markers that have a p-value < return.thresh'<sup>1</sup>
            <h4>log2 fold change</h4>
            The log2 fold change describes a log2 ratio between two groups. In this case the higher this value the better, since it indicates that there is a greater difference between cell clusters 
            regarding the expression of the respective marker gene. 
            The log2 fold change parameter 'limit[s]  testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells'<sup>1</sup>.
            <h4>Download of 3D plots</h4>
            The 3D plots can only be downloaded via the camera icon on the top right corner of every plot.
            <br><br>
            <sup>1</sup> <a href='https://cran.r-project.org/web/packages/Seurat/Seurat.pdf' target='_blank'>Seurat Package Manual</a>
           </div>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  ## Information Button Markergenes
  observeEvent(input$showInfoMark, {
    showModal(modalDialog(
      title=HTML("<b>Additional information:</b>"),
      HTML("<div style='line-height: 1.7'><h4>Limits of the legend in the feature plots</h4>
            - absolute = overall minimum and maximum of the expression of all genes plotted<br>
            - relative = minimum and maximum of the gene expression of every individual gene
            <h4>List of marker genes<sup>1</sup></h4>
            <table style='width:100%;'>
             <tr>
               <td>- p_val</td>
               <td>&nbsp;=&nbsp;</td>
               <td> unadjusted p-value</td>
             </tr>
             <tr>
               <td>- avg_logFC</td>
               <td>&nbsp;=&nbsp;</td>
               <td> log fold-chage of the average expression between the two groups.<br>Positive values indicate that the feature is more highly expressed in the first group.</td>
             </tr>
             <tr>
               <td>- pct.1</td>
               <td>&nbsp;=&nbsp;</td>
               <td> percentage of cells where the feature is detected in the first group</td>
             </tr>
             <tr>
               <td>- pct.2</td>
               <td>&nbsp;=&nbsp;</td>
               <td> percentage of cells where the feature is detected in the second group</td>
             </tr>
             <tr>
               <td>- p_val_adj</td>
               <td>&nbsp;=&nbsp;</td>
               <td> adjusted p-value, based on bonferroni correction using all features in the dataset</td>
             </tr>
            </table>
            <br><br>
            <sup>1</sup> <a href='https://satijalab.org/seurat/v3.0/de_vignette.html' target='_blank'>Seurat - Differential expression testing</a>
           </div>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  ## Information Button Results
  observeEvent(input$showInfoResults, {
    showModal(modalDialog(
      title=HTML("<b>Additional information:</b>"),
      HTML("<div style='line-height: 1.7'>
              <h4>Downloads</h4>
              Every plot can be downloaded in the state it's currently displayed via the first download button.<br>
              <b>Caution:</b> The 3D plots can only be downloaded via the camera icon on the top right corner of every plot.<br>
              Besides that the tables containing the normalization and marker gene values can also be saved on your local computer.<br>
              The following data can be found in the table for the marker genes:<sup>1</sup><br>
              <table style='width:100%;'>
               <tr>
                 <td>- p_val</td>
                 <td>&nbsp;=&nbsp;</td>
                 <td> unadjusted p-value</td>
               </tr>
               <tr>
                 <td>- avg_logFC</td>
                 <td>&nbsp;=&nbsp;</td>
                 <td> log fold-chage of the average expression between the two groups.<br>Positive values indicate that the feature is more highly expressed in the first group.</td>
               </tr>
               <tr>
                 <td>- pct.1</td>
                 <td>&nbsp;=&nbsp;</td>
                 <td> percentage of cells where the feature is detected in the first group</td>
               </tr>
               <tr>
                 <td>- pct.2</td>
                 <td>&nbsp;=&nbsp;</td>
                 <td> percentage of cells where the feature is detected in the second group</td>
               </tr>
               <tr>
                 <td>- p_val_adj</td>
                 <td>&nbsp;=&nbsp;</td>
                 <td> adjusted p-value, based on bonferroni correction using all features in the dataset</td>
               </tr>
              </table>
              At last a summary where all important parameters of the analysis are listed can be downloaded.
              <h5><i>Significant principal components</i></h5>
              'To overcome the extensive technical noise in any single gene for scRNA-seq data, [cells are clustered] based on their PCA scores, 
              with each PC essentially representing a ‘metagene’ that combines information across a correlated gene set.'<sup>2</sup><br>
              'We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.'<sup>3</sup><br>
              This parameter has been calculated with the data shown in the elbow plot.<br>
              <h5><i>Resolution for clustering</i></h5>
              '[The] resolution parameter [...] sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. 
              We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.'<sup>2</sup>
              <h5><i>p-value</i></h5>
              Simply put, in this case the p-value describes the probability of a gene not being a significant marker gene. The lower the p-value the better. 
              The threshold you can set will define which markers will be looked at in the analysis: 'Only return markers that have a p-value < return.thresh'<sup>4</sup>
              <h5><i>log2 fold change</i></h5>
              The log2 fold change describes a log2 ratio between two groups. In this case the higher this value the better, since it indicates that there is a greater difference between cell clusters 
              regarding the expression of the respective marker gene. 
              The log2 fold change parameter 'limit[s]  testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells'<sup>4</sup>.
              <h4>Plots of normalization</h4>
              <h5><i>Elbow plot</i></h5>
              'A ranking of principle components based on the percentage of variance explained by each one.'<sup>2</sup>
              The elbow in this plot indicates the last significant principal component.
              <h5><i>Jackstraw plot</i></h5>
              'A visualization tool for comparing the distribution of p-values for each [principal component] with a uniform distribution (dashed line).'<sup>2</sup>
              <h5><i>VizPCA</i></h5>
              'Visualize top genes associated with reduction components.'<sup>4</sup> (principal components)             
              <h5><i>PC Heatmap</i></h5>
              'Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.'<sup>4</sup><br>
              '[It also] allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses.'<sup>2</sup>
              <h4>Limits of the legend in the feature plots</h4>
              - absolute = overall minimum and maximum of the expression of all genes plotted<br>
              - relative = minimum and maximum of the gene expression of every individual gene
              <br><br>
              <sup>1</sup> <a href='https://satijalab.org/seurat/v3.0/de_vignette.html' target='_blank'>Seurat - Differential expression testing</a> <br>
              <sup>2</sup> <a href='https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html' target='_blank'>Seurat - Guided Clustering Tutorial V3.0</a> <br>
              <sup>3</sup> <a href='https://satijalab.org/seurat/v1.4/pbmc3k_tutorial.html' target='_blank'>Seurat - Guided Clustering Tutorial V1.4</a> <br>
              <sup>4</sup> <a href='https://cran.r-project.org/web/packages/Seurat/Seurat.pdf' target='_blank'>Seurat Package Manual</a>
            </div>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  ## Extraction of csv-file containing expression matrix
  uploadFile <- reactive({
    if(is.null(input$exprmatrix)) {
      return()
    }
	expMatrixFile <- req(input$exprmatrix)
    read.csv(expMatrixFile$datapath, sep=",", row.names = 1)
  })
  
  output$expr_upload_ok <- renderText({
    if(is.null(input$exprmatrix)) {
      return()
    }
    "<i class='glyphicon glyphicon-ok-sign' style='color: #228B22;'></i> <span>Upload of expression matrix successful.</span>"
  })
  
  output$ui.annoDropdown <- renderUI({
    if(is.null(input$annotable)) {
      return()
    }
    annoTable <- read.csv(input$annotable$datapath, sep=",")
    choice <- c(colnames(annoTable))
    selectizeInput("annoParam", "Please select the parameter of the annotation table for the analysis.",
                   choices = choice,
                   selected = choice[2])
  })
  
  output$anno_upload_ok <- renderText({
    if(is.null(input$annotable)) {
      return()
    }
    "<i class='glyphicon glyphicon-ok-sign' style='color: #228B22;'></i> <span>Upload of annotation table successful.</span>"
    })
  
  ## selection of analysis type
  output$ui.radio <- renderUI({
    if(is.null(uploadFile())) {
      return()
    } else {
      return(radioButtons("automan", "Type of Analysis:", choices = list("Automatic", "Manual")))
    }
  })
  
  output$ui.action <- renderUI({
    if(is.null(uploadFile())) {
      return()
    } else {
      return(actionButton("start", "Continue", icon = icon("forward", lib="glyphicon")))
    }
  })
  
  expMatrix <- eventReactive(input$start, {
    uploadFile()
  })
  
  ## creation/extraction of annotation table
  annoTable <- eventReactive(input$start, {
    annotationTable <- input$annotable
    expMatrix <- uploadFile()
    if(is.null(annotationTable)) {
      barcode <- c(colnames(expMatrix))
      celltype_name <- "1"
      celltype <- rep(celltype_name, length(barcode))
      newAnnoTable <- as.data.frame(cbind(barcode, celltype), stringsAsFactors=FALSE)
      parameter <- "celltype"
      return(list(annoTable=newAnnoTable, parameter=parameter))
    } else {
      newAnnoTable <- read.csv(annotationTable$datapath, sep=",")
      parameter <- input$annoParam
      return(list(annoTable=newAnnoTable, parameter=parameter))
    }
  })
  
  
  ### data analysis
  
  ## creating SCE-object and filtering the data 
  build_and_filter <- eventReactive(input$start, {
    exp_matrix <- expMatrix()
    anno <- annoTable()$annoTable
    exp_matrix[is.na(exp_matrix)] <- 0
    
    sce <- SingleCellExperiment(
      assays = list(counts = as.matrix(exp_matrix)),
      colData = anno
    )
      
    expr_genes <- rowSums(counts(sce) > 0) > 0
    sce <- sce[expr_genes, ]
      
	print("Start QC Metrics")
	tmp <- perCellQCMetrics(sce)
	sce$total_counts <- tmp$sum
	sce$total_features_by_counts <- tmp$detected
	rm(tmp)
    print("Finish QC Metrics")
	
    # automatic filtering process
    if(input$automan=="Automatic") {
      
      withProgress(message="Filtering Cells.", detail="This might take some time.", value = 1, {

        # UMIs
		filter_by_total_counts <- (sce$total_counts > (quantile(sce$total_counts, 0.2, na.rm=TRUE)))
        
        # genes
		filter_by_total_features <- (sce$total_features_by_counts > (quantile(sce$total_features_by_counts, 0.2, na.rm=TRUE)))
        
        # valid cells
        sce$valid_cells <- (
          filter_by_total_features &
            filter_by_total_counts
        )
        
        # gene filter
        filter_genes_by_expr <- apply(
          counts(sce[ , colData(sce)$valid_cells]), 
          1, 
          function(x) length(x[x > 1]) >= 2)
        
        rowData(sce)$valid_cells <- filter_genes_by_expr
        dim(sce[rowData(sce)$valid_cells, colData(sce)$valid_cells])
        reducedDim(sce) <- NULL
        
        # merge filtered data (valid cells and genes)
        sce.qc <- sce[rowData(sce)$valid_cells, colData(sce)$valid_cells] # all valid genes & cells (outliers removed)
      })
      
      combined_sces <- list(sce=sce, sce.qc=sce.qc)
      
      return(combined_sces)
      
    } else if(input$automan=="Manual") {
      return(sce)
    }
  })
  
  
  ## manual filtering tab
  output$filter_tab <- renderMenu({
    if(is.null(input$automan) && is.null(input$start)) {
      return()
    } else if(input$automan=="Automatic" && input$start) {
      return()
    } else if(input$automan=="Manual" && input$start) {
      sidebarMenu(menuItem("Filtering", tabName = "filter", icon = icon("filter", lib="glyphicon")))
    }
  })
  
  ## switching to filter tab
  observeEvent(input$start, {
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      return(updateTabItems(session, "sidebar", selected="filter"))
    }
  })
  
  ## UMI slider
  output$ui.sliderUmi <- renderUI({
    sce <- build_and_filter()
	sliderInput("umiCutoff", "UMI Cutoff", min=0, max=max(sce$total_counts), 
                value=quantile(sce$total_counts, 0.01, na.rm=TRUE))			
  })
  
  ## gene slider
  output$ui.sliderGene <- renderUI({
    sce <- build_and_filter()
    sliderInput("geneCutoff", "Gene Cutoff", min=0, max=max(sce$total_features_by_counts), 
                value=quantile(sce$total_features_by_counts, 0.01, na.rm=TRUE))
  })
  
  ## manual cell filtering process
  cellFilter <- reactive({
    sce <- build_and_filter()
    
    if(is.null(input$umiCutoff) && is.null(input$geneCutoff)) {
      return()
    }
    # UMI filter
    filter_by_total_counts <- (sce$total_counts > input$umiCutoff)
    # gene filter
    filter_by_total_features <- (sce$total_features_by_counts > input$geneCutoff)
    
    # valid cells (combination of UMIs and genes)
    sce$valid_cells <- (
      filter_by_total_features &
        filter_by_total_counts
    )
    return(sce)
  })
  
  output$ui.numTranscripts <- renderUI({
    sce <- cellFilter()
    if(is.null(cellFilter())) {
      return()
    }
    numericInput("numTranscripts", label="Minimum transcripts per gene", value=2, min=1, max=max(counts(sce[ , colData(sce)$valid_cells])))
  })
  
  output$ui.numCells <- renderUI({
    sce <- cellFilter()
    if(is.null(cellFilter())) {
      return()
    }
    numericInput("numCells", label="Minimum cells expressing each gene", value=2, min=1, max=sum(sce$valid_cells))
  })
  
  ## manual gene filtering process
  geneFilter <- reactive({
    sce <- cellFilter()
    if(is.null(input$numTranscripts) || is.null(input$numCells)) {
      return()
    }
    filter_genes_by_expr <- apply(
      counts(sce[ , colData(sce)$valid_cells]), 
      1, 
      function(x) length(x[x >= input$numTranscripts]) >= input$numCells)
    
    rowData(sce)$valid_cells <- filter_genes_by_expr
    dim(sce[rowData(sce)$valid_cells, colData(sce)$valid_cells]) 
    reducedDim(sce) <- NULL
    
    # merge filtered data (valid cells and genes)
    sce.qc <- sce[rowData(sce)$valid_cells, colData(sce)$valid_cells]
    return(list(sce=sce, sce.qc=sce.qc))
  })
  
  output$cellsremain <- renderText({
    sce <- cellFilter()
    if(is.null(cellFilter())) {
      return()
    }
    return(paste("Remaining cells after filtering: ", sum(sce$valid_cells), sep=""))
  })
  
  output$genesremain <- renderText({
    sce <- geneFilter()$sce
    if(is.null(geneFilter()$sce)) {
      return()
    }
    return(paste("Remaining genes after filtering: ", sum(rowData(sce)$valid_cells), " out of ", sum(table(rowData(sce)$valid_cells)), sep=""))
  })
  
  output$ui.startNorm <- renderUI({
    if(is.null(geneFilter())) {
      return()
    }
    actionButton("goToNorm", "Filter and start Normalization", icon = icon("forward", lib="glyphicon"))
  })

  
  ## normalization (manual and automatic)
  seuratAnalysis <- reactive({
    sce.qc <- NULL
    if(is.null(input$automan) || is.null(build_and_filter())) {
      return()
    } else if(input$automan=="Automatic" && !is.null(build_and_filter()$sce.qc)) {
      sce.qc <- build_and_filter()$sce.qc
    } else if(input$automan=="Manual" && !is.null(input$goToNorm) && !is.null(geneFilter()$sce.qc)) {
      sce.qc <- geneFilter()$sce.qc
    } else {
      return()
    }
    
	seuratObj <- as.Seurat(sce.qc, data = "counts", project = "Sample")
    
    withProgress(message="Normalizing.", detail="This might take some time.", value = 1, {

      seuratObj <- NormalizeData(
        object = seuratObj,
        normalization.method = "LogNormalize",
        display.progress = FALSE)
	  
	  seuratObj <- FindVariableFeatures(
        object = seuratObj,
		selection.method = "vst",
        mean.function = ExpMean,
        dispersion.function = LogVMR,
        do.plot = FALSE,
        display.progress = FALSE)

	  # feature = all.genes might take very long, in this case omit this parameter, then only the default 2,000 most variable features will be used
	  seuratObj <- ScaleData(
        object = seuratObj,
        display.progress = FALSE)
	  
	  
      seuratObj <- RunPCA(
        object = seuratObj,
        pc.genes = seuratObj@var.genes,
        do.print = FALSE)
      
      seuratObj <- JackStraw(
        seuratObj,
        num.replicate = 100)
		
      seuratObj <- ScoreJackStraw(
        seuratObj,
        dims = 1:20)	  
    })
    return(seuratObj)
  })
  
  ## create normalization matrix
  normMatrix <- reactive({
      seuratObj <- seuratAnalysis()
	  norm_matrix <- as.data.frame(as.matrix(GetAssayData(seuratObj)))
      return(norm_matrix) 
  })
  
  
  ## estimation of significant principal components
  findSignificantPCs <- reactive({
    seuratObj <- seuratAnalysis()
	
	sdev_pca <- Stdev(seuratObj, reduction = "pca")
	
    pc_num <- c(1:length(sdev_pca))
    
    x <- pc_num
    y <- sdev_pca
    pc_df <- data.frame(x, y)
    
    # line connecting endpoints of curve
    p1 <- which.min(x) 
    p2 <- which.max(x) 
    slope <- (y[p2] - y[p1]) / (x[p2] - x[p1]) 
    int <- y[p1] - slope*x[p1] 
    
    # for every point on the curve (xi, pi), the perpendicular line that goes through that point has
    perpslope <- -1/slope 
    perpint <- y - perpslope*x
    
    # the intersection of the perp line(s) with the connecting line is
    xcross <- (int - perpint) / (perpslope - slope) 
    ycross <- slope*xcross + int 
    
    # the distance between the intersection and the point(s) is
    dists <- sqrt((x - xcross)^2 + (y - ycross)^2)
    
    # the index of the farthest point
    elbowp <- which.max(dists)
    
    elbow_params <- list(pc_df=pc_df, p1=p1, p2=p2, elbowp=elbowp, xcross=xcross, ycross=ycross)
    return(elbow_params)
  })
  
  ## switching to normalization tab
  switchToNorm <- eventReactive(input$goToNorm, {
    if(is.null(seuratAnalysis()) || input$automan=="Automatic" || is.null(findSignificantPCs())) {
      return()
    }
    return(updateTabItems(session, "sidebar", selected="normalize"))
  })
  
  ## manual normalization tab
  output$norm_tab <- renderMenu({
    if(is.null(input$goToNorm) || is.null(switchToNorm())) {
      return()
    } else {
      sidebarMenu(menuItem("Normalization", tabName = "normalize", icon = icon("wrench", lib="glyphicon")))
    }
  })
  
  output$ui.significantPCs <- renderUI({
    elbowp <- findSignificantPCs()$elbowp
    pc_num <- findSignificantPCs()$pc_df$x
    if(is.null(findSignificantPCs())) {
      return()
    }
    numericInput("numPCs", label="Number of significant principal components", value=elbowp, min=1, max=max(pc_num))
  })
  
  output$ui.resolution <- renderUI({
	numericInput("numResolution", label="Resolution for clustering", value=0.5, min=0.1, max=2.0, step=0.1)
  })
  
  output$ui.startCluster <- renderUI({
    if(is.null(input$numPCs) || is.null(input$numResolution)) {
      column(7, tags$div(class="norm_button", actionButton("showInfoNorm", "Get Information", icon=icon("info-sign", lib="glyphicon"))))
    } else {
      column(7, tags$div(class="norm_button", actionButton("showInfoNorm", "Get Information", icon=icon("info-sign", lib="glyphicon")),
                         actionButton("goToCluster", "Start Clustering", icon = icon("forward", lib="glyphicon"))))
    }
  })
 
   
  ## plots of normalization
  
  output$elbow_plot <- renderPlot({
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(elbow())
    })
  })
  
  output$jackstraw_plot <- renderPlot({
    if(is.null(seuratAnalysis())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      return(JackStrawPlot(object=seuratAnalysis(), dims = 1:15))
    })
  })
  
  output$vizpca_plot <- renderPlot({
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
	  VizDimLoadings(object = seuratAnalysis(), dims=1:12, nfeatures = 10, ncol=4,  reduction = "pca")
    })
  })
  
  output$pcheatmap <- renderPlot({
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
	  DimHeatmap(object = seuratAnalysis(), dims = 1:15, nfeatures = 10)
    })
  })

  ## clustering (automatic and manual) and calculation of markergenes (automatic)
  manSeuratClustering <- eventReactive(input$goToCluster, {
    seuratObj <- seuratAnalysis()
    elbowp <- input$numPCs
    pc_num <- findSignificantPCs()$pc_df$x
    resol <- input$numResolution
    
    withProgress(message="Clustering Cells.", detail="This might take some time.", value = 1, {

	  seuratObj <- FindNeighbors(
	    seuratObj, 
		dims = 1:pc_num[elbowp])
		
	  seuratObj <- FindClusters(
        seuratObj,
        reduction = "pca",
        dims = 1:pc_num[elbowp],
        resolution = resol,
        verbose = FALSE,
        save.SNN = TRUE)
      
      seuratObj <- RunTSNE(
        object = seuratObj,
		perplexity = (length(colnames(GetAssayData(seuratObj)))/5),
		
        dims = 1:pc_num[elbowp],
        do.fast = TRUE)
      
      seuratObj_3D <- RunTSNE(
        object = seuratObj,
		perplexity = (length(colnames(GetAssayData(seuratObj)))/5),
		
        dims = 1:pc_num[elbowp],
        do.fast = TRUE,
        dim.embed = 3)
      
      seuratObj_UMAP <- RunUMAP(
        object = seuratObj,
		perplexity = (length(colnames(GetAssayData(seuratObj)))/5),
		
        dims = 1:6,
        do.fast = TRUE)
      
      seuratObj_UMAP_3D <- RunUMAP(
        object = seuratObj,
		perplexity = (length(colnames(GetAssayData(seuratObj)))/5),
		
        dims = 1:6,
        do.fast = TRUE,
		n.components = 3L)
    })
    clusteringObj <- list(seuratObj=seuratObj, seuratObj_3D=seuratObj_3D, seuratObj_UMAP=seuratObj_UMAP, seuratObj_UMAP_3D=seuratObj_UMAP_3D)
    return(clusteringObj)
  })
  
  ## switching to cluster tab
  switchToCluster <- eventReactive(input$goToCluster, {
    if(is.null(manSeuratClustering()) || input$automan=="Automatic") {
      return()
    }
    return(updateTabItems(session, "sidebar", selected="cluster"))
  })
  
  ## manual cluster tab
  output$cluster_tab <- renderMenu({
    if(is.null(manSeuratClustering()) || is.null(switchToCluster())) {
      return()
    } else {
      sidebarMenu(menuItem("Clustering", tabName = "cluster", icon = icon("stats", lib="glyphicon")))
    }
  })
  
  ## plots of clustered data
  
  output$tsne_clustered <- renderPlot({
    if(is.null(tsne_clusters())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(tsne_clusters())
    })
  })
  
  output$umap_clustered <- renderPlot({
    if(is.null(umap_clusters())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(umap_clusters())
    })
  })
  
  output$pca_clustered <- renderPlot({
    if(is.null(pca_clusters())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(pca_clusters())
    })
  })
  
  output$mds_clustered <- renderPlot({
    if(is.null(mds_clusters())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(mds_clusters())
    })
  })
  
  output$tsne_clust3d <- renderPlotly({
    if(is.null(tsne_clusters3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      tsne_clusters3d()
    })
  })
  
  output$umap_clust3d <- renderPlotly({
    if(is.null(umap_clusters3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      umap_clusters3d()
    })
  })
  
  output$pca_clust3d <- renderPlotly({
    if(is.null(pca_clusters3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      pca_clusters3d()
    })
  })
  
  output$mds_clust3d <- renderPlotly({
    if(is.null(mds_clusters3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      mds_clusters3d()
    })
  })
  
  
  ## widgets for the calculation of marker genes
  
  output$ui.findMarkers <- renderUI({
    seuratObj <- manSeuratClustering()$seuratObj
	
	clusters <- levels(seuratObj)
    selectInput("find", label = "Find markergenes for the following cluster:", 
                choices = c("None", clusters), 
                selected = clusters[1])
  })
  
  output$ui.compareMarkers <- renderUI({
    seuratObj <- manSeuratClustering()$seuratObj
	
	clusters <- levels(seuratObj)
    find <- input$find
    cluster_names <- NULL
    
    if(is.null(input$find) || is.null(clusters) || is.null(seuratObj)) {
      return()
    }
    
    cluster_names <- unlist(append(cluster_names, sapply(clusters, function(x) if(x!=find) {x})))

    if(input$find=="None") {
      return()
    }
    selectInput("compare", label = "Compare the selected cluster to the following:", 
                choices = c("All", cluster_names), 
                selected = "All")
  })
  
  output$ui.startMarker <- renderUI({
    if(is.null(input$find) || is.null(input$compare)) {
      column(7, tags$div(class="clust_button", actionButton("showInfoClust", "Get Information", icon=icon("info-sign", lib="glyphicon"))))
    } else {
      column(7, tags$div(class="clust_button", actionButton("showInfoClust", "Get Information", icon=icon("info-sign", lib="glyphicon")),
                         actionButton("goToMarker", "Find Markergenes", icon = icon("forward", lib="glyphicon"))))      
    }
  })
  
  
  ## calculation of marker genes (manual)
  seuratMarkers <- eventReactive(input$goToMarker, {
    seuratObj <- manSeuratClustering()$seuratObj
    findMarker <- input$find
    pval <- input$numPval
    logfc <- input$numLogfc
    if(input$compare=="All") {
      compareMarker <- NULL
    } else {
      compareMarker <- input$compare
    }
    
    withProgress(message="Finding Markergenes.", detail="This might take some time.", value = 1, {
      if(input$allMarkers) {
        markers <- FindAllMarkers(
          seuratObj,
          only.pos = TRUE,
          min.pct = 0.25,
          logfc.threshold = logfc,
          return.thresh = pval)
      } else {
        markers <- NULL
      }
      
      if(findMarker!="None") {
        markers1 <- FindMarkers(seuratObj, findMarker, compareMarker, logfc.threshold = logfc)
      } else {
        markers1 <- NULL
      }
    })
    
    if(is.null(markers) && is.null(markers1)) {
      return("Missing input. Could not detect any markers.")
    } else {
      return(list(markers=markers, markers1=markers1))
    }
    
  })

  ## switching to markergene tab
  switchToMarker <- eventReactive(input$goToMarker, {
    if(is.null(seuratMarkers()) || input$automan=="Automatic") {
      return()
    } else {
      updateTabItems(session, "sidebar", selected="marker")
    }
  })
  
  ## manual marker gene tab
  output$marker_tab <- renderMenu({
    if(is.null(input$goToMarker) || is.null(switchToMarker())) {
      return()
    } else {
      sidebarMenu(menuItem("Markergenes", tabName = "marker", icon = icon("tasks", lib="glyphicon")))
    }
  })
  
  ## table of marker genes
  output$ui.downloadMarkersAll <- renderUI({
    if(input$allMarkers==TRUE) {
      downloadButton("downloadAllMarkers", "Download Markergenes for all Clusters")
    }
  })
  
  output$ui.downloadMarkersSpec <- renderUI({
    if(input$find!="None") {
      downloadButton("downloadSpecMarkers", paste("Dowload Markergenes for Cluster", input$find, sep=""))
    }
  })
  
  output$ui.tabsetMarkers <- renderUI({
    if(input$allMarkers && input$find!="None") {
      tabsetPanel(type="tabs",
                  tabPanel("All Clusters",
                           tabsetPanel(type="tabs",
                             tabPanel("Heatmap Top Markers", plotOutput("heatmap1", height="550px")),
                             tabPanel("Expression of Markergenes",
                                      fluidRow(column(6, uiOutput("ui.markers1")), 
                                               column(3, radioButtons("limits", "Limits FeaturePlot Legende",
                                                                      choices = list("Absolute", "Relative"))),
                                               column(3, radioButtons("display", "Display FeaturePlots as...",
                                                                      choices = list("t-SNE", "UMAP")))),
                                      tabsetPanel(type="tabs",
                                                  tabPanel("FeaturePlot", plotOutput("feature1", height="550px")),
                                                  tabPanel("ViolinPlot", plotOutput("violinMarkers1", height="550px")),
                                                  tabPanel("RidgePlot", plotOutput("ridgeMarkers1", height="550px")),
                                                  tabPanel("DotPlot", plotOutput("dot1", height="550px")),
                                                  tabPanel("Heatmap", plotOutput("heatmapMarkers1", height="550px")),
                                                  tabPanel("Co-Expression", plotOutput("coExpr1", height="550px")))
                                      ),
                             tabPanel("List of Markergenes",
                                      fluidRow(column(12, dataTableOutput("tableMarkers1")))))),
                  tabPanel(paste("Cluster", input$find, sep=""),
                           tabsetPanel(type="tabs",
                                       tabPanel("Heatmap Top Markers", plotOutput("heatmap2", height="550px")),
                                       tabPanel("Expression of Markergenes", 
                                                fluidRow(column(6, uiOutput("ui.markers2")),
                                                         column(3, radioButtons("limits2", "Limits FeaturePlot Legende",
                                                                                choices = list("Absolute", "Relative"))),
                                                         column(3, radioButtons("display2", "Display FeaturePlots as...",
                                                                                choices = list("t-SNE", "UMAP")))),
                                                tabsetPanel(type="tabs",
                                                            tabPanel("FeaturePlot", plotOutput("feature2", height="550px")),
                                                            tabPanel("ViolinPlot", plotOutput("violinMarkers2", height="550px")),
                                                            tabPanel("RidgePlot", plotOutput("ridgeMarkers2", height="550px")),
                                                            tabPanel("DotPlot", plotOutput("dot2", height="550px")),
                                                            tabPanel("Heatmap", plotOutput("heatmapMarkers2", height="550px")),
                                                            tabPanel("Co-Expression", plotOutput("coExpr2", height="550px")))
                                                ),
                                       tabPanel("List of Markergenes",
                                                fluidRow(column(12, dataTableOutput("tableMarkers2"))))))
      )
    } else if(input$allMarkers && input$find=="None") {
      tabsetPanel(type="tabs",
                  tabPanel("Heatmap all Clusters", plotOutput("heatmap1", height="550px")),
                  tabPanel("Expression of Markergenes all Clusters",
                           fluidRow(column(6, uiOutput("ui.markers1")),
                                    column(3, radioButtons("limits", "Limits FeaturePlot Legende",
                                                           choices = list("Absolute", "Relative"))),
                                    column(3, radioButtons("display", "Display FeaturePlots as...",
                                                           choices = list("t-SNE", "UMAP")))),
                           tabsetPanel(type="tabs",
                                       tabPanel("FeaturePlot", plotOutput("feature1", height="550px")),
                                       tabPanel("ViolinPlot", plotOutput("violinMarkers1", height="550px")),
                                       tabPanel("RidgePlot", plotOutput("ridgeMarkers1", height="550px")),
                                       tabPanel("DotPlot", plotOutput("dot1", height="550px")),
                                       tabPanel("Heatmap", plotOutput("heatmapMarkers1", height="550px")),
                                       tabPanel("Co-Expression", plotOutput("coExpr1", height="550px")))
                           ),
                  tabPanel("List of Markergenes all Clusters",
                           fluidRow(column(12, dataTableOutput("tableMarkers1")))))
    } else if(input$allMarkers==FALSE && input$find!="None") {
      tabsetPanel(type="tabs",
                  tabPanel(paste("Heatmap Top Markers Cluster", input$find, sep=""), plotOutput("heatmap2", height="550px")),
                  tabPanel(paste("Expression of Markergenes Cluster", input$find, sep=""),
                           fluidRow(column(6, uiOutput("ui.markers2")),
                                    column(3, radioButtons("limits2", "Limits FeaturePlot Legende",
                                                           choices = list("Absolute", "Relative"))),
                                    column(3, radioButtons("display2", "Display FeaturePlots as...",
                                                           choices = list("t-SNE", "UMAP")))),
                           tabsetPanel(type="tabs",
                                       tabPanel("FeaturePlot", plotOutput("feature2", height="550px")),
                                       tabPanel("ViolinPlot", plotOutput("violinMarkers2", height="550px")),
                                       tabPanel("RidgePlot", plotOutput("ridgeMarkers2", height="550px")),
                                       tabPanel("DotPlot", plotOutput("dot2", height="550px")),
                                       tabPanel("Heatmap", plotOutput("heatmapMarkers2", height="550px")),
                                       tabPanel("Co-Expression", plotOutput("coExpr2", height="550px")))
                           ),
                  tabPanel(paste("List of Markergenes Cluster", input$find, sep=""),
                           fluidRow(column(12, dataTableOutput("tableMarkers2")))))
    } else if(input$allMarkers==FALSE && input$find=="None") {
      return(as.character(seuratMarkers()))
    }
  })
  
  ### plots marker genes
  
  ## heatmap for all clusters
  output$heatmap1 <- renderPlot({
    if(is.null(top5_heatmap())) {
      return()
    }
    withProgress(message="Clustering Cells.", detail="This might take some time.", value = 1, {
      plot(top5_heatmap())
    })
  })
  
  ## heatmap for specific cluster
  output$heatmap2 <- renderPlot({
    if(is.null(top20_heatmap())) {
      return()
    }
    withProgress(message="Clustering Cells.", detail="This might take some time.", value = 1, {
      plot(top20_heatmap())
    })
  })

  ## dropdown menu for all clusters
  output$ui.markers1 <- renderUI({
    markers <- seuratMarkers()$markers
    selectizeInput("markergenes1", "Markergenes (max. 3):",
                   multiple = TRUE,
                   options = list(maxItems = 3),
                   choices= c((as.character(markers[order(markers$avg_logFC, decreasing = TRUE),]$gene))),
                   selected = c(markers[order(markers$avg_logFC, decreasing = TRUE),]$gene[1:3]))
  })

  ## dropdown menu for specific cluster
  output$ui.markers2 <- renderUI({
    markers <- seuratMarkers()$markers1
    markers$gene <- rownames(markers)
    selectizeInput("markergenes2", "Markergenes (max. 3):",
                   multiple = TRUE,
                   options = list(maxItems = 3),
                   choices= c((as.character(markers[order(markers$avg_logFC, decreasing = TRUE),]$gene))),
                   selected = c(markers[order(markers$avg_logFC, decreasing = TRUE),]$gene[1:3]))
  })

  ## feature plot for all clusters
  output$feature1 <- renderPlot({
    seuratObj <- manSeuratClustering()$seuratObj
    seuratObj_UMAP <- manSeuratClustering()$seuratObj_UMAP
    markers <- input$markergenes1
    
    if(input$limits=="Absolute") {
	  limits <- c(min(GetAssayData(object = seuratObj)), max(GetAssayData(object = seuratObj)))
    } else {
      limits <- c()
    }
    
    if(input$display=="t-SNE") {
      reduction <- "tsne"
      name <- "t-SNE"
    } else {
      reduction <- "umap"
      name <- "UMAP"
      seuratObj <- seuratObj_UMAP
    }
    
    if(is.null(markers)) {
      return()
    } else {
      withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
        fp <- FeaturePlot(
          object = seuratObj,
		  features = markers,
		  cols=c("lightgrey", "#3c4cbd"),
          ncol=length(markers),
		  label.size = 4,
		  pt.size = 2.5,
		  reduction = reduction,
		  # Important new parameter as Seurat 3 stores FeaturePlots differently per default, causing the lapply modification of each single plot to crash.
		  combine = FALSE
        )
		
        fps <- lapply(fp, function(x) {x + labs(x=paste(name, "1", sep="_"), y=paste(name, "2", sep="_")) + 
            scale_colour_gradientn(limits = limits, colours=c("#2e6c8f", "#3c8dbc", "#bc3c4d", "#8f2e3b"), name="gene exp") +
            theme(plot.title=element_text(face="plain", size=18), 
                  axis.text=element_text(size=12), 
                  axis.title=element_text(face="plain", size=14),
                  legend.title=element_text(face="plain", size=12))})

				  
        print(cowplot::plot_grid(plotlist = fps, ncol = length(markers)))
      })
    }
  })
  
  ## feature plot for specific cluster
  output$feature2 <- renderPlot({
    seuratObj <- manSeuratClustering()$seuratObj
    seuratObj_UMAP <- manSeuratClustering()$seuratObj_UMAP
    markers <- input$markergenes2
    
    if(input$limits2=="Absolute") {
	  limits <- c(min(GetAssayData(object = seuratObj)), max(GetAssayData(object = seuratObj)))
    } else {
      limits <- c()
    }
    
    if(input$display2=="t-SNE") {
      reduction <- "tsne"
      name <- "t-SNE"
    } else {
      reduction <- "umap"
      name <- "UMAP"
      seuratObj <- seuratObj_UMAP
    }
    
    if(is.null(markers)) {
      return()
    } else {
      withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
		fp <- FeaturePlot(
          object = seuratObj,
		  features = markers,
		  cols=c("lightgrey", "#3c4cbd"),
          ncol=length(markers),
		  label.size = 4,
		  pt.size = 2.5,
		  reduction = reduction,
		  combine = FALSE
        )
        fps <- lapply(fp, function(x) {x + labs(x=paste(name, "1", sep="_"), y=paste(name, "2", sep="_")) + 
            scale_colour_gradientn(limits = limits, colours=c("#2e6c8f", "#3c8dbc", "#bc3c4d", "#8f2e3b"), name="gene exp") +
            theme(plot.title=element_text(face="plain", size=18), 
                  axis.text=element_text(size=12), 
                  axis.title=element_text(face="plain", size=14),
                  legend.title=element_text(face="plain", size=12))
        })
        cowplot::plot_grid(plotlist = fps, ncol = length(markers))
      })
    }
  })

  
  ## Violin Plot Function (customizing Violin plot)
  
  # Function to change x labels of individual violin plots
  VlnPlot_2 <- function(object, features.plot, xlab, nCol, colours, qc) {
    # Main function
    main_function <- function(object = object, features.plot = features.plot, xlab = xlab, colour = colours) {
      VlnPlot(object = object, features = features.plot, cols = colour ) + #size.x.use = 14, size.y.use = 14, size.title.use = 18
        labs(x = xlab, y = "Frequency") + 
        theme(plot.title = element_text(face = "plain", size=18), axis.text = element_text(size=12), 
              axis.title.x =element_text(face="plain", colour="black", size=14), axis.title.y = element_text(face = "plain", colour = "black", size=14))
    }
	
	qc_function <- function(object = object, features.plot = features.plot, xlab = xlab, colour = colours) {
      VlnPlot(object = object, features = features.plot, cols = colour ) + #size.x.use = 14, size.y.use = 14, size.title.use = 18
        labs(x = xlab, y = "Frequency") + 
        theme(plot.title = element_text(face = "plain", size=18), axis.text = element_text(size=12), 
              axis.title.x =element_text(face="plain", colour="black", size=14), axis.title.y = element_text(face = "plain", colour = "black", size=14), legend.position="none") + {
			  if(features.plot == "total_features_by_counts") ggtitle("Total Features") else ggtitle("Total Counts")
			  }
    }
    
    # Apply main function on all features
	if(qc == FALSE){
		p <- lapply(X = features.plot, object = object, xlab = xlab, FUN = main_function)
	} else if(qc == TRUE) {
		p <- lapply(X = features.plot, object = object, xlab = xlab, FUN = qc_function)
	}
    
    # Arrange all plots using cowplot
    cowplot::plot_grid(plotlist = p, ncol = nCol)
  }
  
  ## violin plot for all clusters
  violin1 <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      seuratObj <- seuratClustering()$seuratObj
      markers <- input$markergenes
    } else {
      seuratObj <- manSeuratClustering()$seuratObj
      markers <- input$markergenes1
    }
    
    if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    VlnPlot_2(object = seuratObj, features.plot = c(markers), xlab = "Cluster", nCol = 1, colours = colour, qc = FALSE)
  })
  
  output$violinMarkers1 <- renderPlot({
    seuratObj <- manSeuratClustering()$seuratObj
    markers <- input$markergenes1
    
    if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      VlnPlot_2(object = seuratObj, features.plot = c(markers), xlab = "Cluster", nCol = length(markers), colours = colour, qc = FALSE)
    })
  })

  ## violin plot for specific cluster
  violin2 <- reactive({
    seuratObj <- manSeuratClustering()$seuratObj
    markers <- input$markergenes2
    
    if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    VlnPlot_2(object = seuratObj, features.plot = c(markers), xlab = "Cluster", nCol = 1, colours = colour, qc = FALSE)
  })
  
  output$violinMarkers2 <- renderPlot({
    seuratObj <- manSeuratClustering()$seuratObj
    markers <- input$markergenes2
    
    if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      VlnPlot_2(object = seuratObj, features.plot = c(markers), xlab = "Cluster", nCol = length(markers), colours = colour, qc = FALSE)
    })
  })
  
  ## ridge plot function
  RidgePlot_2 <- function(object, features.plot, nCol, colours) {
    
    # Main function
    main_function <- function(object = object, features.plot = features.plot, colour = colours) {
      RidgePlot(object = object, features = features.plot, cols = colour) +
        labs(x = "Frequency", y = "Cluster") + 
        theme(plot.title = element_text(face = "plain", hjust = 0.5, size=18), axis.text = element_text(size=12), 
              axis.title.x =element_text(face="plain", colour="black", hjust = 0.5, size=14), axis.title.y = element_text(face = "plain", colour = "black", hjust = 0.5, size=14))
    }
	
    # Apply main function on all features
    p <- lapply(X = features.plot, object = object, FUN = main_function)
    
    # Arrange all plots using cowplot
    cowplot::plot_grid(plotlist = p, ncol = nCol)
  }
  
  ## ridge plot for all clusters
  ridge1 <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      seuratObj <- seuratClustering()$seuratObj
      markers <- input$markergenes
    } else {
      seuratObj <- manSeuratClustering()$seuratObj
      markers <- input$markergenes1
    }
    
    if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    RidgePlot_2(object = seuratObj, features.plot = c(markers), nCol = 1, colours = colour) 
  })
  
  output$ridgeMarkers1 <- renderPlot({
    seuratObj <- manSeuratClustering()$seuratObj
    markers <- input$markergenes1
    
    if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      RidgePlot_2(object = seuratObj, features.plot = c(markers), nCol = length(markers), colours = colour)
    })
  })
  
  ## ridge plot for specific cluster
  ridge2 <- reactive({
    seuratObj <- manSeuratClustering()$seuratObj
    markers <- input$markergenes2
    
    if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    RidgePlot_2(object = seuratObj, features.plot = c(markers), nCol = 1, colours = colour)
  })
  
  output$ridgeMarkers2 <- renderPlot({
    seuratObj <- manSeuratClustering()$seuratObj
    markers <- input$markergenes2
    
    if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      RidgePlot_2(object = seuratObj, features.plot = c(markers), nCol = length(markers), colours = colour)
    })
  })
  
  ## dot plot for all clusters
  dotplot1 <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      seuratObj <- seuratClustering()$seuratObj
      markers <- input$markergenes
    } else {
      seuratObj <- manSeuratClustering()$seuratObj
      markers <- input$markergenes1
    }
    
    if(is.null(markers)) {
      return()
    }

	DotPlot(object = seuratObj, features = markers, cols = c("#2e6c8f", "#8f2e3b")) + ggtitle("Gene expression per cluster") +
      scale_colour_gradientn(colours=c("#2e6c8f", "#3c8dbc", "#bc3c4d", "#8f2e3b")) + labs(size="% cells", color="avg expr") +
      theme(plot.title=element_text(face="plain", size=18), 
            axis.text=element_text(size=12), 
            axis.title=element_text(face="plain", size=14),
            legend.title=element_text(face="plain", size=12))
    
  })
  
  output$dot1 <- renderPlot({
    if(is.null(dotplot1())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(dotplot1())
    })
  })
  
  ## dot plot for specific cluster
  dotplot2 <- reactive({
    seuratObj <- manSeuratClustering()$seuratObj
    markers <- input$markergenes2

    if(is.null(markers)) {
      return()
    }

	DotPlot(object = seuratObj, features = markers, cols = c("#2e6c8f", "#8f2e3b")) + ggtitle("Gene expression per cluster") +
        scale_colour_gradientn(colours=c("#2e6c8f", "#3c8dbc", "#bc3c4d", "#8f2e3b")) + labs(size="% cells", color="avg expr") +
        theme(plot.title=element_text(face="plain", size=18), 
                axis.text=element_text(size=12), 
                axis.title=element_text(face="plain", size=14),
                legend.title=element_text(face="plain", size=12))
  })
  
  output$dot2 <- renderPlot({
    if(is.null(dotplot2())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(dotplot2())
    })
  })
  
  ## heatmap with selected marker genes for all clusters
  heatMark1 <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      seuratObj <- seuratClustering()$seuratObj
      markers <- input$markergenes
    } else {
      seuratObj <- manSeuratClustering()$seuratObj
      markers <- input$markergenes1
    }
    
    if(is.null(markers)) {
      return()
    }
	
    DoHeatmap(
      object = seuratObj,
      features = markers)	
  })
  
  output$heatmapMarkers1 <- renderPlot({
    if(is.null(heatMark1())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(heatMark1())
    })
  })
  
  ## heatmap with selected marker genes for specific cluster
  heatMark2 <- reactive({
    seuratObj <- manSeuratClustering()$seuratObj
    markers <- input$markergenes2
    
    if(is.null(markers)) {
      return()
    }

    DoHeatmap(
      object = seuratObj,
      features = markers)
  })
  
  output$heatmapMarkers2 <- renderPlot({
    if(is.null(heatMark2())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(heatMark2())
    })
  })
  
  ## Co-Exression FeaturePlot
  
  ## function for customizing the feature plot (co-expression) of Seurat
  CoExprPlot_2 <- function(fpBoth, fpLegend, sObj, name, geneX, geneY) {
    main_function_both <- function(fp=fp, sObj1=sObj, name1=name) {
      fp <- fp + labs(x=paste(name1, "1", sep="_"), y=paste(name1, "2", sep="_")) + theme(plot.title=element_text(face="plain", size=18), 
                                                                                          axis.text=element_text(size=12), 
                                                                                          axis.title=element_text(face="plain", size=14),
																						  legend.position="none") #legend.title=element_text(face="plain", size=12),
    }
	
	fp_both <- fpBoth + labs(x=paste(name, "1", sep="_"), y=paste(name, "2", sep="_")) + theme(plot.title=element_text(face="plain", size=18), 
                                                                                          axis.text=element_text(size=12), 
                                                                                          axis.title=element_text(face="plain", size=14),
																					      legend.position="none")
	fp_legend <- fpLegend + labs(x=geneX, y=geneY) + theme(plot.title=element_text(face="plain", size=18), 
                                                                                          axis.text=element_text(size=12), 
                                                                                          axis.title=element_text(face="plain", size=14),
																						  legend.position="none")
	
	cowplot::plot_grid(fp_both, fp_legend, ncol = 2)
    
  }
  
  ## feature plot for co-expression of marker genes for all clusters
  coEx1 <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      seuratObj <- seuratClustering()$seuratObj
      seuratObj_UMAP <- seuratClustering()$seuratObj_UMAP
      markers <- input$markergenes[1:2]
    } else {
      seuratObj <- manSeuratClustering()$seuratObj
      seuratObj_UMAP <- manSeuratClustering()$seuratObj_UMAP
      markers <- input$markergenes1[1:2]
    }
    
    if(input$display=="t-SNE") {
      reduction <- "tsne"
      name <- "t-SNE"
    } else {
      reduction <- "umap"
      name <- "UMAP"
      seuratObj <- seuratObj_UMAP
    }
    
    if(!is.na(markers[2]) && !is.null(markers)) {
	  fp <- FeaturePlot(object = seuratObj, features = markers, cols=c("lightgrey", "#2e3b8f", "#3b8f2e", "#8f2e3b"), 
            reduction = reduction, pt.size=2.5, blend = TRUE)
      fp <- CoExprPlot_2(fp[[3]], fp[[4]], seuratObj, name, markers[1], markers[2])
      return(fp)
    } else {
      return(NULL)
    }
    
  })
  
  output$coExpr1 <- renderPlot({
    if(is.null(coEx1())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      coEx1()
    })
  })
  
  ## feature plot for co-expression of marker genes for specific cluster
  coEx2 <- reactive({
    seuratObj <- manSeuratClustering()$seuratObj
    seuratObj_UMAP <- manSeuratClustering()$seuratObj_UMAP
    markers <- input$markergenes2[1:2]

    if(input$display=="t-SNE") {
      reduction <- "tsne"
      name <- "t-SNE"
    } else {
      reduction <- "umap"
      name <- "UMAP"
      seuratObj <- seuratObj_UMAP
    }
    
    if(!is.na(markers[2]) && !is.null(markers)) {
							
	  fp <- FeaturePlot(object = seuratObj, features = markers, cols=c("lightgrey", "#2e3b8f", "#3b8f2e", "#8f2e3b"), 
            reduction = reduction, pt.size=2.5, blend = TRUE)
    
	  fp <- CoExprPlot_2(fp[[3]], fp[[4]], seuratObj, name, markers[1], markers[2])
      return(fp)
    } else {
      return(NULL)
    }
  })
  
  output$coExpr2 <- renderPlot({
    if(is.null(coEx2())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      coEx2()
    })
  })
  
  ## table of marker genes (all clusters)
  output$tableMarkers1 <- DT::renderDataTable(
    DT::datatable(seuratMarkers()$markers, rownames = TRUE, options = list(pageLength = 20))
  )
  
  ## table of marker genes (specific cluster)
  output$tableMarkers2 <- DT::renderDataTable(
    DT::datatable(seuratMarkers()$markers1, rownames=TRUE, options = list(pageLength = 20))
  )
  
  
  ### clustering and calculation of marker gens (automatic analysis)
  
  seuratClustering <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic" && !is.null(seuratAnalysis()) && !is.null(findSignificantPCs())) {
      seuratObj <- seuratAnalysis()
      elbowp <- findSignificantPCs()$elbowp
      pc_num <- findSignificantPCs()$pc_df$x
      
      withProgress(message="Clustering Cells.", detail="This might take some time.", value = 1, {

	  seuratObj <- FindNeighbors(
	    seuratObj, 
		dims = 1:pc_num[elbowp])

        seuratObj <- FindClusters(
          seuratObj,
          reduction.type = "pca",
          dims.use = 1:pc_num[elbowp],
		  resolution = 0.5,
          verbose = FALSE,
          save.SNN = TRUE)
		
        seuratObj <- RunTSNE(
          object = seuratObj,
          perplexity = (length(colnames(GetAssayData(seuratObj)))/5),
          dims = 1:pc_num[elbowp],
          do.fast = TRUE)
        
        seuratObj_3D <- RunTSNE(
          object = seuratObj,
          perplexity = (length(colnames(GetAssayData(seuratObj)))/5),
          dims = 1:pc_num[elbowp],
          do.fast = TRUE,
          dim.embed = 3)
        
        seuratObj_UMAP <- RunUMAP(
          object = seuratObj,
          perplexity = (length(colnames(GetAssayData(seuratObj)))/5),
          dims = 1:6,
          do.fast = TRUE)
        
        seuratObj_UMAP_3D <- RunUMAP(
          object = seuratObj,
          perplexity = (length(colnames(GetAssayData(seuratObj)))/5),
          dims = 1:6,
          do.fast = TRUE,
          n.components = 3L)
        
        markers <- FindAllMarkers(
          seuratObj,
          only.pos = TRUE,
          min.pct = 0.25,
          return.thresh = 0.01)
      })
      clusteringObj <- list(seuratObj=seuratObj, seuratObj_3D=seuratObj_3D, seuratObj_UMAP=seuratObj_UMAP, seuratObj_UMAP_3D=seuratObj_UMAP_3D, markers=markers)
      return(clusteringObj)
    } else {
      return()
    }
  })

  
  ### results
  
  output$ui.downloadButtons <- renderUI({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      tags$div(downloadButton("report", "Download All Plots"),
      downloadButton("normtable", "Download Normalization Table"),
      downloadButton("markerlist", "Download List of Markergenes"),
      downloadButton("logfile", "Download Summary"))
    } else {
      if(input$allMarkers && input$find!="None") {
        tags$div(downloadButton("report", "Download All Plots"),
                 downloadButton("normtable", "Download Normalization Table"),
                 downloadButton("downloadAllMarkers1", "Download Markergenes for all Clusters"),
                 downloadButton("downloadSpecMarkers1", paste("Dowload Markergenes for Cluster", input$find, sep="")),
                 downloadButton("logfile", "Download Summary"))
      } else if(input$allMarkers && input$find=="None") {
        tags$div(downloadButton("report", "Download All Plots"),
                 downloadButton("normtable", "Download Normalization Table"),
                 downloadButton("downloadAllMarkers1", "Download Markergenes for all Clusters"),
                 downloadButton("logfile", "Download Summary"))
      } else if(input$allMarkers==FALSE && input$find!="None") {
        tags$div(downloadButton("report", "Download All Plots"),
                 downloadButton("normtable", "Download Normalization Table"),
                 downloadButton("downloadSpecMarkers1", paste("Dowload Markergenes for Cluster", input$find, sep="")),
                 downloadButton("logfile", "Download Summary"))
      } else if(input$allMarkers==FALSE && input$find=="None") {
        tags$div(downloadButton("report", "Download All Plots"),
                 downloadButton("normtable", "Download Normalization Table"),
                 downloadButton("logfile", "Download Summary"))
      }
    }
  })
  
  output$ui.tabsetResults <- renderUI({
    if(is.null(input$automan)) {
      return() 
    } else if(input$automan=="Manual" && input$allMarkers && input$find!="None") {
      tabsetPanel(type="tabs",
                  tabPanel("Quality Control",
                           tabsetPanel(type="tabs",
                                       tabPanel("Filtered by Counts", plotOutput("plot1", height="550px")),
                                       tabPanel("Filtered by Features", plotOutput("plot2", height="550px")),
                                       tabPanel("Gene Distribution", plotOutput("plot2_2", height="550px")),
                                       tabPanel("Violin Plot", plotOutput("plot6", height="550px")))),
                  tabPanel("Raw Data", 
                           tabsetPanel(type="tabs",
                                       tabPanel("2D",
                                                checkboxInput("scaleCells_raw", label = "Scale cells by number of genes", value = TRUE),
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotOutput("plot3", height="550px")),
                                                            tabPanel("UMAP", plotOutput("plot3_1", height="550px")),
                                                            tabPanel("PCA", plotOutput("plot3_2", height="550px")),
                                                            tabPanel("MDS", plotOutput("plot3_3", height="550px")))),
                                       tabPanel("3D",
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotlyOutput("plot3_4", height="550px")),
                                                            tabPanel("UMAP", plotlyOutput("plot3_5", height="550px")),
                                                            tabPanel("PCA", plotlyOutput("plot3_6", height="550px")),
                                                            tabPanel("MDS", plotlyOutput("plot3_7", height="550px"))))
                                       )
                           ),
                  tabPanel("Filtered Data", 
                           tabsetPanel(type="tabs",
                                       tabPanel("2D",
                                                checkboxInput("scaleCells_valid", label = "Scale cells by number of genes", value = TRUE),
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotOutput("plot4", height="550px")),
                                                            tabPanel("UMAP", plotOutput("plot4_1", height="550px")),
                                                            tabPanel("PCA", plotOutput("plot4_2", height="550px")),
                                                            tabPanel("MDS", plotOutput("plot4_3", height="550px")))),
                                       tabPanel("3D",
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotlyOutput("plot4_4", height="550px")),
                                                            tabPanel("UMAP", plotlyOutput("plot4_5", height="550px")),
                                                            tabPanel("PCA", plotlyOutput("plot4_6", height="550px")),
                                                            tabPanel("MDS", plotlyOutput("plot4_7", height="550px"))))
                                       )
                           ),
                  tabPanel("Normalized Data", 
                           tabsetPanel(type="tabs",
                                       tabPanel("2D",
                                                checkboxInput("scaleCells_norm", label = "Scale cells by number of genes", value = TRUE),
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotOutput("plot5", height="550px")),
                                                            tabPanel("UMAP", plotOutput("plot5_1", height="550px")),
                                                            tabPanel("PCA", plotOutput("plot5_2", height="550px")),
                                                            tabPanel("MDS", plotOutput("plot5_3", height="550px")))),
                                       tabPanel("3D",
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotlyOutput("plot5_4", height="550px")),
                                                            tabPanel("UMAP", plotlyOutput("plot5_5", height="550px")),
                                                            tabPanel("PCA", plotlyOutput("plot5_6", height="550px")),
                                                            tabPanel("MDS", plotlyOutput("plot5_7", height="550px"))))
                                       )
                           ),
                  tabPanel("Normalization",
                           tabsetPanel(type="tabs",
                                       tabPanel("Elbow Plot", plotOutput("plot7", height="550px")),
                                       tabPanel("Jackstraw Plot", plotOutput("plot8", height="550px")),
                                       tabPanel("VizPCA Plot", plotOutput("plot8_2", height="550px")),
                                       tabPanel("PCHeatmap", plotOutput("plot8_3", height="550px")))),
                  tabPanel("Clustered Cells",
                           tabsetPanel(type="tabs",
                                       tabPanel("2D",
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotOutput("plot9", height="550px")),
                                                            tabPanel("UMAP", plotOutput("plot9_1", height="550px")),
                                                            tabPanel("PCA", plotOutput("plot9_2", height="550px")),
                                                            tabPanel("MDS", plotOutput("plot9_3", height="550px")))),
                                       tabPanel("3D",
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotlyOutput("plot9_4", height="550px")),
                                                            tabPanel("UMAP", plotlyOutput("plot9_5", height="550px")),
                                                            tabPanel("PCA", plotlyOutput("plot9_6", height="550px")),
                                                            tabPanel("MDS", plotlyOutput("plot9_7", height="550px"))))
                                       )
                           ),
                  tabPanel("All Clusters",
                           tabsetPanel(type="tabs",
                                       tabPanel("Heatmap Top Markers", plotOutput("plot10", height="550px")),
                                       tabPanel("FeaturePlot", plotOutput("plot11", height="550px")),
                                       tabPanel("ViolinPlot", plotOutput("plot12", height="550px")),
                                       tabPanel("RidgePlot", plotOutput("plot12_2", height="550px")),
                                       tabPanel("DotPlot", plotOutput("plot12_3", height="550px")),
                                       tabPanel("Heatmap", plotOutput("plot10_2", height="550px")),
                                       tabPanel("Co-Expression", plotOutput("plot11_2", height="550px"))
                           )),
                  tabPanel(paste("Cluster", input$find, sep=""),
                           tabsetPanel(type="tabs",
                                       tabPanel("Heatmap Top Markers", plotOutput("plot13", height="550px")),
                                       tabPanel("FeaturePlot", plotOutput("plot14", height="550px")),
                                       tabPanel("ViolinPlot", plotOutput("plot15", height="550px")),
                                       tabPanel("RidgePlot", plotOutput("plot15_2", height="550px")),
                                       tabPanel("DotPlot", plotOutput("plot15_3", height="550px")),
                                       tabPanel("Heatmap", plotOutput("plot13_2", height="550px")),
                                       tabPanel("Co-Expression", plotOutput("plot14_2", height="550px"))
                           ))
      )
    } else {
      tabsetPanel(type="tabs",
                  tabPanel("Quality Control",
                           tabsetPanel(type="tabs",
                                       tabPanel("Filtered by Counts", plotOutput("plot1", height="550px")),
                                       tabPanel("Filtered by Features", plotOutput("plot2", height="550px")),
                                       tabPanel("Gene Distribution", plotOutput("plot2_2", height="550px")),
                                       tabPanel("Violin Plot", plotOutput("plot6", height="550px")))),
                  tabPanel("Raw Data", 
                           tabsetPanel(type="tabs",
                                       tabPanel("2D",
                                                checkboxInput("scaleCells_raw", label = "Scale cells by number of genes", value = TRUE),
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotOutput("plot3", height="550px")),
                                                            tabPanel("UMAP", plotOutput("plot3_1", height="550px")),
                                                            tabPanel("PCA", plotOutput("plot3_2", height="550px")),
                                                            tabPanel("MDS", plotOutput("plot3_3", height="550px")))),
                                       tabPanel("3D",
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotlyOutput("plot3_4", height="550px")),
                                                            tabPanel("UMAP", plotlyOutput("plot3_5", height="550px")),
                                                            tabPanel("PCA", plotlyOutput("plot3_6", height="550px")),
                                                            tabPanel("MDS", plotlyOutput("plot3_7", height="550px"))))
                                       )
                           ),
                  tabPanel("Filtered Data", 
                           tabsetPanel(type="tabs",
                                       tabPanel("2D",
                                                checkboxInput("scaleCells_valid", label = "Scale cells by number of genes", value = TRUE),
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotOutput("plot4", height="550px")),
                                                            tabPanel("UMAP", plotOutput("plot4_1", height="550px")),
                                                            tabPanel("PCA", plotOutput("plot4_2", height="550px")),
                                                            tabPanel("MDS", plotOutput("plot4_3", height="550px")))),
                                       tabPanel("3D",
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotlyOutput("plot4_4", height="550px")),
                                                            tabPanel("UMAP", plotlyOutput("plot4_5", height="550px")),
                                                            tabPanel("PCA", plotlyOutput("plot4_6", height="550px")),
                                                            tabPanel("MDS", plotlyOutput("plot4_7", height="550px"))))
                                       )
                           ),
                  tabPanel("Normalized Data", 
                           tabsetPanel(type="tabs",
                                       tabPanel("2D",
                                                checkboxInput("scaleCells_norm", label = "Scale cells by number of genes", value = TRUE),
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotOutput("plot5", height="550px")),
                                                            tabPanel("UMAP", plotOutput("plot5_1", height="550px")),
                                                            tabPanel("PCA", plotOutput("plot5_2", height="550px")),
                                                            tabPanel("MDS", plotOutput("plot5_3", height="550px")))),
                                       tabPanel("3D",
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotlyOutput("plot5_4", height="550px")),
                                                            tabPanel("UMAP", plotlyOutput("plot5_5", height="550px")),
                                                            tabPanel("PCA", plotlyOutput("plot5_6", height="550px")),
                                                            tabPanel("MDS", plotlyOutput("plot5_7", height="550px"))))
                                       )
                           ),
                  tabPanel("Normalization",
                           tabsetPanel(type="tabs",
                                       tabPanel("Elbow Plot", plotOutput("plot7", height="550px")),
                                       tabPanel("Jackstraw Plot", plotOutput("plot8", height="550px")),
                                       tabPanel("VizPCA Plot", plotOutput("plot8_2", height="550px")),
                                       tabPanel("PCHeatmap", plotOutput("plot8_3", height="550px")))),
                  tabPanel("Clustered Cells",
                           tabsetPanel(type="tabs",
                                       tabPanel("2D",
                                                if(input$automan=="Automatic") {
                                                  checkboxInput("scaleCells_clust", label = "Scale cells by number of genes", value = TRUE)
                                                },
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotOutput("plot9", height="550px")),
                                                            tabPanel("UMAP", plotOutput("plot9_1", height="550px")),
                                                            tabPanel("PCA", plotOutput("plot9_2", height="550px")),
                                                            tabPanel("MDS", plotOutput("plot9_3", height="550px")))),
                                       tabPanel("3D",
                                                tabsetPanel(type="tabs",
                                                            tabPanel("t-SNE", plotlyOutput("plot9_4", height="550px")),
                                                            tabPanel("UMAP", plotlyOutput("plot9_5", height="550px")),
                                                            tabPanel("PCA", plotlyOutput("plot9_6", height="550px")),
                                                            tabPanel("MDS", plotlyOutput("plot9_7", height="550px"))))
                                       )
                           ),
                  if(input$automan=="Automatic") {
                    tabPanel("Expression of Markergenes",
                             fluidRow(
                               column(6, uiOutput("ui.markers")),
                               column(3, radioButtons("limits", "Limits FeaturePlot Legend",
                                                      choices = list("Absolute", "Relative"))),
                               column(3, radioButtons("display", "Display FeaturePlots as...",
                                                      choices = list("t-SNE", "UMAP")))
                             ),
                             tabsetPanel(type="tabs",
                                         tabPanel("FeaturePlot", plotOutput("plot11", height="550px")),
                                         tabPanel("ViolinPlot", plotOutput("plot12", height="550px")),
                                         tabPanel("RidgePlot", plotOutput("plot12_2", height="550px")),
                                         tabPanel("DotPlot", plotOutput("plot12_3", height="550px")),
                                         tabPanel("Heatmap", plotOutput("plot10_2", height="550px")),
                                         tabPanel("Co-Expression", plotOutput("plot11_2", height="550px")),
                                         tabPanel("Heatmap Top Markergenes", plotOutput("plot10", height="550px")))
                    )
                  } else if(input$automan=="Manual" && input$allMarkers && input$find=="None") {
                    tabPanel("All Clusters",
                             tabsetPanel(type="tabs",
                                         tabPanel("Heatmap Top Markers", plotOutput("plot10", height="550px")),
                                         tabPanel("FeaturePlot", plotOutput("plot11", height="550px")),
                                         tabPanel("ViolinPlot", plotOutput("plot12", height="550px")),
                                         tabPanel("RidgePlot", plotOutput("plot12_2", height="550px")),
                                         tabPanel("DotPlot", plotOutput("plot12_3", height="550px")),
                                         tabPanel("Heatmap", plotOutput("plot10_2", height="550px")),
                                         tabPanel("Co-Expression", plotOutput("plot11_2", height="550px"))
                             ))
                  } else if(input$automan=="Manual" && input$allMarkers==FALSE && input$find!="None") {
                    tabPanel(paste("Cluster", input$find, sep=""),
                             tabsetPanel(type="tabs",
                                         tabPanel("Heatmap Top Markers", plotOutput("plot13", height="550px")),
                                         tabPanel("FeaturePlot", plotOutput("plot14", height="550px")),
                                         tabPanel("ViolinPlot", plotOutput("plot15", height="550px")),
                                         tabPanel("RidgePlot", plotOutput("plot15_2", height="550px")),
                                         tabPanel("DotPlot", plotOutput("plot15_3", height="550px")),
                                         tabPanel("Heatmap", plotOutput("plot13_2", height="550px")),
                                         tabPanel("Co-Expression", plotOutput("plot14_2", height="550px"))
                             ))
                  } else if(input$automan=="Manual"&& input$allMarkers==FALSE && input$find=="None") {
                    tabPanel("Markergenes", as.character(seuratMarkers()))
                  }
      )
    }
  })
  
  ## switching to results tab (manual analysis)
  switchToResult <- eventReactive(input$goToResult, {
    updateTabItems(session, "sidebar", selected="results")
  })
  
  ## result tab for automatic and manual analysis
  output$result_tab <- renderMenu({
    if(is.null(input$automan)) {
      return()
    } else if((!is.null(seuratClustering()) && input$automan=="Automatic") || (!is.null(switchToResult()) && input$automan=="Manual")){
      sidebarMenu(menuItem("Results", tabName = "results", icon = icon("download-alt", lib="glyphicon"))) 
    } else {
      return()
    }
  })

  ## switching to results tab (automatic analysis)
  observe({
    if(is.null(seuratClustering()) || input$automan!="Automatic") {
      return()
    } else {
      return(updateTabItems(session, "sidebar", selected="results"))
    }
  })

  ### plots for every step in the analysis
  
  ## cutoff in UMI histogram
  abline_tc <- reactive({
    if(input$automan=="Automatic") {
      sce <- build_and_filter()$sce
      quant_tc <- quantile(sce$total_counts, 0.05, na.rm=TRUE)
      return(quant_tc)
    } else if(input$automan=="Manual") {
      sce <- build_and_filter()
      cutoffUmi <- input$umiCutoff
      return(cutoffUmi)
    }
  })
  
  ## cutoff in gene histogram
  abline_tf <- reactive({
    if(input$automan=="Automatic") {
      sce <- build_and_filter()$sce
      quant_tf <- quantile(sce$total_features_by_counts, 0.05, na.rm=TRUE)
      return(quant_tf)
    } else if(input$automan=="Manual") {
      sce <- build_and_filter()
      cutoffGene <- input$geneCutoff
      return(cutoffGene)
    }
  })
  
  ## UMI histogram
  histo_tc <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce <- build_and_filter()
    } else {
      sce <- build_and_filter()$sce
    }
    line_tc <- abline_tc()
    histo <- ggplot(data=as.data.frame(sce$total_counts), aes(sce$total_counts)) + 
      geom_histogram(bins = length(sce$total_counts),
                     fill=I("white"),
                     col=I("black")) + 
      geom_vline(xintercept = as.numeric(line_tc), colour="red") +
      ggtitle("Histogram of UMI counts per cell") +
      labs(x="UMI counts", y="cell counts") +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, size=18), 
            axis.text = element_text(size=12), axis.title = element_text(size=14))
  })
  
  output$umihisto <- renderPlot({
    if(is.null(histo_tc())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(histo_tc())
    })
  })

  ## gene histogram
  histo_tf <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce <- build_and_filter()
    } else {
      sce <- build_and_filter()$sce
    }
    line_tf <- abline_tf()
    ggplot(data=as.data.frame(sce$total_features_by_counts), aes(sce$total_features_by_counts)) + 
      geom_histogram(bins = length(sce$total_features_by_counts),
                     fill=I("white"),
                     col=I("black")) + 
      geom_vline(xintercept = as.numeric(line_tf), colour="red") +
      ggtitle("Histogram of gene counts per cell") +
      labs(x="gene counts", y="cell counts") +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, size=18), 
            axis.text = element_text(size=12), axis.title = element_text(size=14))
  })
  
  output$genehisto <- renderPlot({
    if(is.null(histo_tf())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(histo_tf())
    })
  })
  
  ## scatter plot for gene filter
  scatter_genes <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce <- build_and_filter()
    } else {
      sce <- build_and_filter()$sce
    }
    
    mean_exp_per_gene <- rowMeans(counts(sce)) # average expression of genes
    cells_per_gene <- rowSums(counts(sce) > 0) # number of cells expressing the specific gene
    df_genefilter <- cbind(mean_exp_per_gene, cells_per_gene)
    
    ggplot(data=as.data.frame(df_genefilter), aes(x=mean_exp_per_gene, y=cells_per_gene)) + 
      geom_point() + geom_smooth(se = FALSE, color="red") +
      theme_classic(base_size=14) + ggtitle("Frequency distribution of genes") + labs(x="Mean transcripts per gene", y="Frequency of genes") +
      theme(plot.title = element_text(hjust = 0.5, size=18), 
            axis.text = element_text(size=12), axis.title = element_text(size=14))
  })
  
  output$genedistrib <- renderPlot({
    if(is.null(scatter_genes())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(scatter_genes())
    })
  })
  
  ## plots of raw data
  tsne_raw <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce <- geneFilter()$sce
    } else {
      sce <- build_and_filter()$sce
    }
    if(input$scaleCells_raw == TRUE)
      size <- "total_features_by_counts"
    else
      size <- NULL
	sce <- runTSNE(sce, exprs_values = "counts", rand_seed = 123456)
    raw_tsne <- plotTSNE(sce,
                         size_by = size,
                         colour_by = "valid_cells",
                         shape_by = annoTable()$parameter
    ) + ggtitle("t-SNE of raw data") + theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
                                             axis.text = element_text(size=12), axis.title = element_text(size=14))
    raw_tsne$guides$size$title <- "gene counts"
    
    if(length(unique(annoTable()$annoTable[,annoTable()$parameter]))==1) {
      raw_tsne <- raw_tsne + labs(x="t-SNE_1", y="t-SNE_2", fill="valid cells") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc"))
      return(raw_tsne)
    } else {
      raw_tsne <- raw_tsne + labs(x="t-SNE_1", y="t-SNE_2", colour="valid cells") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc"), aesthetics = "colour")
      return(raw_tsne)
    }
  })
  
  umap_raw <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce <- geneFilter()$sce
    } else {
      sce <- build_and_filter()$sce
    }
    if(input$scaleCells_raw == TRUE) {
      size <- "total_features_by_counts"
    } else {
      size <- NULL
    }
	sce <- runUMAP(sce, exprs_values = "counts")
    raw_umap <- plotUMAP(sce,
                         size_by = size,
                         colour_by = "valid_cells",
                         shape_by = annoTable()$parameter
    ) + ggtitle("UMAP of raw data") + theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
                                            axis.text = element_text(size=12), axis.title = element_text(size=14))
    raw_umap$guides$size$title <- "gene counts"
    
    if(length(unique(annoTable()$annoTable[,annoTable()$parameter]))==1) {
      raw_umap <- raw_umap + labs(x="UMAP_1", y="UMAP_2", fill="valid cells") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc"))
      return(raw_umap)
    } else {
      raw_umap <- raw_umap + labs(x="UMAP_1", y="UMAP_2", colour="valid cells") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc"), aesthetics = "colour")
      return(raw_umap)
    }
  })
  
  pca_raw <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce <- geneFilter()$sce
    } else {
      sce <- build_and_filter()$sce
    }
    if(input$scaleCells_raw == TRUE) {
      size <- "total_features_by_counts"
    } else {
      size <- NULL
    }
	sce <- runPCA(sce, exprs_values = "counts")	
    raw_pca <- plotPCASCE(sce,
                          size_by = size,
                          colour_by = "valid_cells",
                          shape_by = annoTable()$parameter
    ) + ggtitle("PCA of raw data") + theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
                                           axis.text = element_text(size=12), axis.title = element_text(size=14))
    raw_pca$guides$size$title <- "gene counts"
    
    if(length(unique(annoTable()$annoTable[,annoTable()$parameter]))==1) {
      raw_pca <- raw_pca + labs(x="PCA_1", y="PCA_2", fill="valid cells") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc"))
      return(raw_pca)
    } else {
      raw_pca <- raw_pca + labs(x="PCA_1", y="PCA_2", colour="valid cells") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc"), aesthetics = "colour")
      return(raw_pca)
    }
  })
  
  mds_raw <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce <- geneFilter()$sce
    } else {
      sce <- build_and_filter()$sce
    }
    
    sce <- runMDS(sce, exprs_values = "counts")
    
    if(input$scaleCells_raw == TRUE) {
      size <- "total_features_by_counts"
    } else {
      size <- NULL
    }
	sce <- runMDS(sce, exprs_values = "counts")
	
    raw_mds <- plotMDS(sce,
                       size_by = size,
                       colour_by = "valid_cells",
                       shape_by = annoTable()$parameter
    ) + ggtitle("MDS of raw data") + theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
                                           axis.text = element_text(size=12), axis.title = element_text(size=14))
    raw_mds$guides$size$title <- "gene counts"
    
    if(length(unique(annoTable()$annoTable[,annoTable()$parameter]))==1) {
      raw_mds <- raw_mds + labs(x="MDS_1", y="MDS_2", fill="valid cells") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc"))
      return(raw_mds)
    } else {
      raw_mds <- raw_mds + labs(x="MDS_1", y="MDS_2", colour="valid cells") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc"), aesthetics = "colour")
      return(raw_mds)
    }
  })
  
  ## 3d plots of raw data
  tsne_raw3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce <- geneFilter()$sce
    } else {
      sce <- build_and_filter()$sce
    }
    
    sce <- runTSNE(sce, ncomponents = 3, exprs_values = "counts", rand_seed = 123456)
	tSNE_tmp <- reducedDim(sce, "TSNE", withDimnames = TRUE)
	colnames(tSNE_tmp) <- c("tSNE_1", "tSNE_2", "tSNE_3")
    rownames(tSNE_tmp) <- c(colnames(counts(sce)))
	reducedDim(sce, "TSNE", withDimnames = TRUE) <- tSNE_tmp
	rm(tSNE_tmp)
    
	plot_ly(data.frame(reducedDim(sce, "TSNE", withDimnames = TRUE)), x=~tSNE_1,y=~tSNE_2, z=~tSNE_3, type="scatter3d", mode="markers", color=~sce$valid_cells, colors=c("#bc6b3c", "#3c8dbc"), 
            size=~sce$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce)), '<br>Gene count:', sce$total_features_by_counts, '<br>UMI count:', sce$total_counts)) %>%
      layout(sce, title="3D t-SNE of raw data")
  })
  
  umap_raw3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce <- geneFilter()$sce
    } else {
      sce <- build_and_filter()$sce
    }
    
    sce <- runUMAP(sce, ncomponents = 3, exprs_values = "counts")
    UMAP_tmp <- reducedDim(sce, "UMAP", withDimnames = TRUE)
	colnames(UMAP_tmp) <- c("UMAP_1", "UMAP_2", "UMAP_3")
	reducedDim(sce, "UMAP", withDimnames = TRUE) <- UMAP_tmp
	rm(UMAP_tmp)
    
    plot_ly(data.frame(reducedDim(sce, "UMAP", withDimnames = TRUE)), x=~UMAP_1,y=~UMAP_2, z=~UMAP_3, type="scatter3d", mode="markers", color=~sce$valid_cells, colors=c("#bc6b3c", "#3c8dbc"), 
            size=~sce$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce)), '<br>Gene count:', sce$total_features_by_counts, '<br>UMI count:', sce$total_counts)) %>%
      layout(sce, title="3D UMAP of raw data")
  })
  
  pca_raw3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce <- geneFilter()$sce
    } else {
      sce <- build_and_filter()$sce
    }
    
	### SCATER UPDATE ###
    sce <- runPCA(sce, ncomponents = 3, exprs_values = "counts")
    PCA_tmp <- reducedDim(sce, "PCA", withDimnames = TRUE)
	colnames(PCA_tmp) <- c("PC_1", "PC_2", "PC_3")
	reducedDim(sce, "PCA", withDimnames = TRUE) <- PCA_tmp
	rm(PCA_tmp)
    
    plot_ly(data.frame(reducedDim(sce, "PCA", withDimnames = TRUE)), x=~PC_1,y=~PC_2, z=~PC_3, type="scatter3d", mode="markers", color=~sce$valid_cells, colors=c("#bc6b3c", "#3c8dbc"), 
            size=~sce$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce)), '<br>Gene count:', sce$total_features_by_counts, '<br>UMI count:', sce$total_counts)) %>%
      layout(sce, title="3D PCA of raw data")
  })
  
  mds_raw3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce <- geneFilter()$sce
    } else {
      sce <- build_and_filter()$sce
    }
    
	### SCATER UPDATE ###
    sce <- runMDS(sce, ncomponents = 3, exprs_values = "counts")
    MDS_tmp <- reducedDim(sce, "MDS", withDimnames = TRUE)
	colnames(MDS_tmp) <- c("MDS_1", "MDS_2", "MDS_3")
	reducedDim(sce, "MDS", withDimnames = TRUE) <- MDS_tmp
	rm(MDS_tmp)
    
    plot_ly(data.frame(reducedDim(sce, "MDS", withDimnames = TRUE)), x=~MDS_1,y=~MDS_2, z=~MDS_3, type="scatter3d", mode="markers", color=~sce$valid_cells, colors=c("#bc6b3c", "#3c8dbc"), 
            size=~sce$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce)), '<br>Gene count:', sce$total_features_by_counts, '<br>UMI count:', sce$total_counts)) %>%
      layout(sce, title="3D MDS of raw data")
  })
  
  ## plots of filtered data
  tsne_validcells <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce.qc <- geneFilter()$sce.qc
    } else {
      sce.qc <- build_and_filter()$sce.qc
    }
    if(input$scaleCells_valid == TRUE) {
      size <- "total_features_by_counts"
    } else {
      size <- NULL
    }
	sce.qc <- runTSNE(sce.qc, exprs_values = "counts", rand_seed = 123456) 
    valid_tsne <- plotTSNE(sce.qc,
                           size_by = size,
                           shape_by = annoTable()$parameter  
    ) + ggtitle("t-SNE of valid cells") + labs(x="t-SNE_1", y="t-SNE_2") + scale_fill_manual(values="#3c8dbc") +
      theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
            axis.text = element_text(size=12), axis.title = element_text(size=14))
    valid_tsne$guides$size$title <- "gene counts"
    
    return(valid_tsne)
  })
  
  umap_validcells <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce.qc <- geneFilter()$sce.qc
    } else {
      sce.qc <- build_and_filter()$sce.qc
    }
    if(input$scaleCells_valid == TRUE) {
      size <- "total_features_by_counts"
    } else {
      size <- NULL
    }
    sce.qc <- runUMAP(sce.qc, exprs_values = "counts")
	valid_umap <- plotUMAP(sce.qc,
                           size_by = size,
                           shape_by = annoTable()$parameter
    ) + ggtitle("UMAP of valid cells") + labs(x="UMAP_1", y="UMAP_2") + scale_fill_manual(values="#3c8dbc") +
      theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
            axis.text = element_text(size=12), axis.title = element_text(size=14))
    valid_umap$guides$size$title <- "gene counts"
    
    return(valid_umap)
  })
  
  pca_validcells <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce.qc <- geneFilter()$sce.qc
    } else {
      sce.qc <- build_and_filter()$sce.qc
    }
    if(input$scaleCells_valid == TRUE) {
      size <- "total_features_by_counts"
    } else {
      size <- NULL
    }
	sce.qc <- runPCA(sce.qc, exprs_values = "counts")
    valid_pca <- plotPCASCE(sce.qc,
                            size_by = size,
                            shape_by = annoTable()$parameter
    ) + ggtitle("PCA of valid cells") + labs(x="PCA_1", y="PCA_2") + scale_fill_manual(values="#3c8dbc") +
      theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
            axis.text = element_text(size=12), axis.title = element_text(size=14))
    valid_pca$guides$size$title <- "gene counts"
    
    return(valid_pca)
  })
  
  mds_validcells <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce.qc <- geneFilter()$sce.qc
    } else {
      sce.qc <- build_and_filter()$sce.qc
    }
    
    sce.qc <- runMDS(sce.qc, exprs_values = "counts")
    
    if(input$scaleCells_valid == TRUE) {
      size <- "total_features_by_counts"
    } else {
      size <- NULL
    }
	sce.qc <- runMDS(sce.qc, exprs_values = "counts")
    valid_mds <- plotMDS(sce.qc,
                         size_by = size,
                         shape_by = annoTable()$parameter
    ) + ggtitle("MDS of valid cells") + labs(x="MDS_1", y="MDS_2") + scale_fill_manual(values="#3c8dbc") +
      theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
            axis.text = element_text(size=12), axis.title = element_text(size=14))
    valid_mds$guides$size$title <- "gene counts"
    
    return(valid_mds)
  })
  
  ## 3d plots of filtered data
  tsne_validcells3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce.qc <- geneFilter()$sce.qc
    } else {
      sce.qc <- build_and_filter()$sce.qc
    }
    
	### SCATER UPDATE ###
    sce.qc <- runTSNE(sce.qc, ncomponents = 3, exprs_values = "counts", rand_seed = 123456)
    tSNE_tmp <- reducedDim(sce.qc, "TSNE", withDimnames = TRUE)
	colnames(tSNE_tmp) <- c("tSNE_1", "tSNE_2", "tSNE_3")
    rownames(tSNE_tmp) <- c(colnames(counts(sce.qc)))
	reducedDim(sce.qc, "TSNE", withDimnames = TRUE) <- tSNE_tmp
	rm(tSNE_tmp)
    
	plot_ly(data.frame(reducedDim(sce.qc, "TSNE", withDimnames = TRUE)), x=~tSNE_1,y=~tSNE_2, z=~tSNE_3, type="scatter3d", mode="markers", color=sce.qc$valid_cells, colors=c("#3c8dbc"),
            size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
      layout(sce.qc, title="3D t-SNE of valid cells")
  })
  
  umap_validcells3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce.qc <- geneFilter()$sce.qc
    } else {
      sce.qc <- build_and_filter()$sce.qc
    }
    
    sce.qc <- runUMAP(sce.qc, ncomponents = 3, exprs_values = "counts")
    UMAP_tmp <- reducedDim(sce.qc, "UMAP", withDimnames = TRUE)
	colnames(UMAP_tmp) <- c("UMAP_1", "UMAP_2", "UMAP_3")
	reducedDim(sce.qc, "UMAP", withDimnames = TRUE) <- UMAP_tmp
	rm(UMAP_tmp)
    
	plot_ly(data.frame(reducedDim(sce.qc, "UMAP", withDimnames = TRUE)), x=~UMAP_1,y=~UMAP_2, z=~UMAP_3, type="scatter3d", mode="markers", color=sce.qc$valid_cells, colors=c("#3c8dbc"),
            size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
      layout(sce.qc, title="3D UMAP of valid cells")
  })
  
  pca_validcells3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce.qc <- geneFilter()$sce.qc
    } else {
      sce.qc <- build_and_filter()$sce.qc
    }
    
    sce.qc <- runPCA(sce.qc, ncomponents = 3, exprs_values = "counts")
    PCA_tmp <- reducedDim(sce.qc, "PCA", withDimnames = TRUE)
	colnames(PCA_tmp) <- c("PC_1", "PC_2", "PC_3")
	reducedDim(sce.qc, "PCA", withDimnames = TRUE) <- PCA_tmp
	rm(PCA_tmp)
    
	plot_ly(data.frame(reducedDim(sce.qc, "PCA", withDimnames = TRUE)), x=~PC_1,y=~PC_2, z=~PC_3, type="scatter3d", mode="markers", color=sce.qc$valid_cells, colors=c("#3c8dbc"),
            size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
      layout(sce.qc, title="3D PCA of valid cells")
  })
  
  mds_validcells3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      sce.qc <- geneFilter()$sce.qc
    } else {
      sce.qc <- build_and_filter()$sce.qc
    }
    
    sce.qc <- runMDS(sce.qc, ncomponents = 3, exprs_values = "counts")
    MDS_tmp <- reducedDim(sce.qc, "MDS", withDimnames = TRUE)
	colnames(MDS_tmp) <- c("MDS_1", "MDS_2", "MDS_3")
	reducedDim(sce.qc, "MDS", withDimnames = TRUE) <- MDS_tmp
	rm(MDS_tmp)
	
	plot_ly(data.frame(reducedDim(sce.qc, "MDS", withDimnames = TRUE)), x=~MDS_1,y=~MDS_2, z=~MDS_3, type="scatter3d", mode="markers", color=sce.qc$valid_cells, colors=c("#3c8dbc"),
            size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
      layout(sce.qc, title="3D MDS of valid cells")
  })
  
  ## plots of normalized data
  tsne_norm <- reactive({
    seuratObj <- seuratAnalysis()
	sce.qc <- as.SingleCellExperiment(seuratObj)
    if(input$scaleCells_norm == TRUE) {
      size <- "total_features_by_counts"
    } else {
      size <- NULL
    }
	sce.qc <- runTSNE(sce.qc, exprs_values = "counts", rand_seed = 123456)
    norm_tsne <- plotTSNE(sce.qc,
                          size_by = size,
                          shape_by = annoTable()$parameter 
    ) + ggtitle("t-SNE of normalized data") + labs(x="t-SNE_1", y="t-SNE_2") + scale_fill_manual(values="#3c8dbc") +
      theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
            axis.text = element_text(size=12), axis.title = element_text(size=14))
    norm_tsne$guides$size$title <- "gene counts"
    
    return(norm_tsne)
  })
  
  umap_norm <- reactive({
    seuratObj <- seuratAnalysis()
	
	sce.qc <- as.SingleCellExperiment(seuratObj)
    if(input$scaleCells_norm == TRUE) {
      size <- "total_features_by_counts"
    } else {
      size <- NULL
    }
	sce.qc <- runUMAP(sce.qc, exprs_values = "counts")
    norm_umap <- plotUMAP(sce.qc,
                          size_by = size,
                          shape_by = annoTable()$parameter
    ) + ggtitle("UMAP of normalized data") + labs(x="UMAP_1", y="UMAP_2") + scale_fill_manual(values="#3c8dbc") +
      theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
            axis.text = element_text(size=12), axis.title = element_text(size=14))
    norm_umap$guides$size$title <- "gene counts"
    
    return(norm_umap)
  })
  
  pca_norm <- reactive({
    seuratObj <- seuratAnalysis()
	sce.qc <- as.SingleCellExperiment(seuratObj)
    if(input$scaleCells_norm == TRUE) {
      size <- "total_features_by_counts"
    } else {
      size <- NULL
    }
	sce.qc <- runPCA(sce.qc, exprs_values = "counts")
    norm_pca <- plotPCASCE(sce.qc,
                           size_by = size,
                           shape_by = annoTable()$parameter
    ) + ggtitle("PCA of normalized data") + labs(x="PCA_1", y="PCA_2") + scale_fill_manual(values="#3c8dbc") +
      theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
            axis.text = element_text(size=12), axis.title = element_text(size=14))
    norm_pca$guides$size$title <- "gene counts"
    
    return(norm_pca)
  })
  
  mds_norm <- reactive({
    seuratObj <- seuratAnalysis()
	sce.qc <- as.SingleCellExperiment(seuratObj)
    sce.qc <- runMDS(sce.qc, exprs_values = "logcounts")
    if(input$scaleCells_norm == TRUE) {
      size <- "total_features_by_counts"
    } else {
      size <- NULL
    }
	sce.qc <- runMDS(sce.qc, exprs_values = "counts")
    norm_mds <- plotMDS(sce.qc,
                        size_by = size,
                        shape_by = annoTable()$parameter
    ) + ggtitle("MDS of normalized data") + labs(x="MDS_1", y="MDS_2") + scale_fill_manual(values="#3c8dbc") +
      theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
            axis.text = element_text(size=12), axis.title = element_text(size=14))
    norm_mds$guides$size$title <- "gene counts"
    
    return(norm_mds)
  })
  
  ## 3d plots of normalized data
  tsne_norm3d <- reactive({
    seuratObj <- seuratAnalysis()
	sce.qc <- as.SingleCellExperiment(seuratObj)
    
    sce.qc <- runTSNE(sce.qc, ncomponents = 3, exprs_values = "logcounts", rand_seed = 123456)
    tSNE_tmp <- reducedDim(sce.qc, "TSNE", withDimnames = TRUE)
	colnames(tSNE_tmp) <- c("tSNE_1", "tSNE_2", "tSNE_3")
    rownames(tSNE_tmp) <- c(colnames(logcounts(sce.qc)))
	reducedDim(sce.qc, "TSNE", withDimnames = TRUE) <- tSNE_tmp
	rm(tSNE_tmp)

	plot_ly(data.frame(reducedDim(sce.qc, "TSNE", withDimnames = TRUE)), x=~tSNE_1,y=~tSNE_2, z=~tSNE_3, type="scatter3d", mode="markers", color=sce.qc$valid_cells, colors=c("#3c8dbc"),
            size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
      layout(sce.qc, title="3D t-SNE of normalized data")
  })
  
  umap_norm3d <- reactive({
    seuratObj <- seuratAnalysis()
	sce.qc <- as.SingleCellExperiment(seuratObj)
    
    sce.qc <- runUMAP(sce.qc, ncomponents = 3, exprs_values = "logcounts")
	UMAP_tmp <- reducedDim(sce.qc, "UMAP", withDimnames = TRUE)
	colnames(UMAP_tmp) <- c("UMAP_1", "UMAP_2", "UMAP_3")
	reducedDim(sce.qc, "UMAP", withDimnames = TRUE) <- UMAP_tmp
	rm(UMAP_tmp)
	
	plot_ly(data.frame(reducedDim(sce.qc, "UMAP", withDimnames = TRUE)), x=~UMAP_1,y=~UMAP_2, z=~UMAP_3, type="scatter3d", mode="markers", color=sce.qc$valid_cells, colors=c("#3c8dbc"),
            size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
      layout(sce.qc, title="3D UMAP of normalized data")
  })
  
  pca_norm3d <- reactive({
    seuratObj <- seuratAnalysis()
	sce.qc <- as.SingleCellExperiment(seuratObj)
    
    sce.qc <- runPCA(sce.qc, ncomponents = 3, exprs_values = "logcounts")
    PCA_tmp <- reducedDim(sce.qc, "PCA", withDimnames = TRUE)
	colnames(PCA_tmp) <- c("PC_1", "PC_2", "PC_3")
	reducedDim(sce.qc, "PCA", withDimnames = TRUE) <- PCA_tmp
	rm(PCA_tmp)
	
	plot_ly(data.frame(reducedDim(sce.qc, "PCA", withDimnames = TRUE)), x=~PC_1,y=~PC_2, z=~PC_3, type="scatter3d", mode="markers", color=sce.qc$valid_cells, colors=c("#3c8dbc"),
            size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
      layout(sce.qc, title="3D PCA of normalized data")
  })
  
  mds_norm3d <- reactive({
    seuratObj <- seuratAnalysis()
	sce.qc <- as.SingleCellExperiment(seuratObj)
    
    sce.qc <- runMDS(sce.qc, ncomponents = 3, exprs_values = "logcounts")
    MDS_tmp <- reducedDim(sce.qc, "MDS", withDimnames = TRUE)
	colnames(MDS_tmp) <- c("MDS_1", "MDS_2", "MDS_3")
	reducedDim(sce.qc, "MDS", withDimnames = TRUE) <- MDS_tmp
	rm(MDS_tmp)
    
	plot_ly(data.frame(reducedDim(sce.qc, "MDS", withDimnames = TRUE)), x=~MDS_1,y=~MDS_2, z=~MDS_3, type="scatter3d", mode="markers", color=sce.qc$valid_cells, colors=c("#3c8dbc"),
            size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
      layout(sce.qc, title="3D MDS of normalized data")
  })
  
  ## plots of normalization
  violin <- reactive({
    seuratObj <- seuratAnalysis()
	VlnPlot_2(object = seuratObj, features.plot = c("total_features_by_counts", "total_counts"), xlab = "", nCol = 2, colours = "#3c8dbc", qc = TRUE)
  })
  
  elbow <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      pc_df <- findSignificantPCs()$pc_df
      elbowp <- findSignificantPCs()$elbowp
    } else if(input$automan=="Manual") {
      if(is.null(seuratAnalysis()) || is.null(findSignificantPCs())) {
        return()
      }
      pc_df <- findSignificantPCs()$pc_df
      elbowp <- findSignificantPCs()$elbowp
    }

    ggplot(data=pc_df, aes(x=pc_df$x, y=pc_df$y)) +
      geom_point(data=pc_df, aes(x=pc_df$x, y=pc_df$y), size=3.0) +
      geom_point(data=pc_df, aes(x=pc_df$x[elbowp], y=pc_df$y[elbowp]), col="red", size=3.0) +
      ggtitle("Elbow Plot of principal components") +
      labs(x="principal component", y="Standard deviation of PC") +
      theme_classic(base_size = 15) +
      theme(plot.title = element_text(hjust = 0.5))
  })

  ## plots of clustered data
  tsne_clusters <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      seuratObj <- manSeuratClustering()$seuratObj
      if(input$scaleCells_clust1 == TRUE) {
        size <- "total_features_by_counts"
      } else {
        size <- NULL
      }
    } else {
      seuratObj <- seuratClustering()$seuratObj
      if(input$scaleCells_clust == TRUE) {
        size <- "total_features_by_counts"
      } else {
        size <- NULL
      }
    }
	sce.qc <- as.SingleCellExperiment(seuratObj)
    sce.qc$clusters <- sce.qc$ident # "clusters" instead of "ident" in plot
    
    tsne <- plotTSNE(sce.qc, 
                     colour_by = "clusters", 
                     size_by = size,
                     shape_by = annoTable()$parameter
    ) + ggtitle("t-SNE of clustered cells") + theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
                                                    axis.text = element_text(size=12), axis.title = element_text(size=14))
    tsne$guides$size$title <- "gene counts"
    
    if(length(levels(sce.qc$clusters)) <= 7 && length(unique(annoTable()$annoTable[,annoTable()$parameter]))==1) {
      tsne <- tsne + labs(x="t-SNE_1", y="t-SNE_2", fill="cluster")  + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc"))
      return(tsne)
    } else if(length(levels(sce.qc$clusters)) <= 7 && length(unique(annoTable()$annoTable[,annoTable()$parameter]))>1) {
      tsne <- tsne + labs(x="t-SNE_1", y="t-SNE_2", colour="cluster") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc"), aesthetics = "colour")
      return(tsne)
    } else {
      return(tsne)
    }
  })
  
  umap_clusters <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      seuratObj <- manSeuratClustering()$seuratObj_UMAP
      if(input$scaleCells_clust1 == TRUE) {
        size <- "total_features_by_counts"
      } else {
        size <- NULL
      }
    } else {
      seuratObj <- seuratClustering()$seuratObj_UMAP
      if(input$scaleCells_clust == TRUE) {
        size <- "total_features_by_counts"
      } else {
        size <- NULL
      }
    }
	sce.qc <- as.SingleCellExperiment(seuratObj)
    sce.qc$clusters <- sce.qc$ident # "clusters" instead of "ident" in plot
    
    umap <-  plotUMAP(object=sce.qc,
                      colour_by = "clusters",
                      size_by = size,
                      shape_by = annoTable()$parameter
    ) + ggtitle("UMAP of clustered cells") + theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"),
                                                   axis.text = element_text(size=12), axis.title = element_text(size=14))
    umap$guides$size$title <- "gene counts"
    
    if(length(levels(sce.qc$clusters)) <= 7 && length(unique(annoTable()$annoTable[,annoTable()$parameter]))==1) {
      umap <- umap + labs(x="UMAP_1", y="UMAP_2", fill="cluster") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")) 
      return(umap)
    } else if(length(levels(sce.qc$clusters)) <= 7 && length(unique(annoTable()$annoTable[,annoTable()$parameter]))>1) {
      umap <- umap + labs(x="UMAP_1", y="UMAP_2", colour="cluster") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc"), aesthetics = "colour")
      return(umap)
    } else {
      return(umap)
    }
  })
  
  pca_clusters <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      seuratObj <- manSeuratClustering()$seuratObj
      if(input$scaleCells_clust1 == TRUE) {
        size <- "total_features_by_counts"
      } else {
        size <- NULL
      }
    } else {
      seuratObj <- seuratClustering()$seuratObj
      if(input$scaleCells_clust == TRUE) {
        size <- "total_features_by_counts"
      } else {
        size <- NULL
      }
    }
	sce.qc <- as.SingleCellExperiment(seuratObj)
    sce.qc$clusters <- sce.qc$ident # "clusters" instead of "ident" in plot
    
    pca <- plotPCASCE(sce.qc, 
                      colour_by = "clusters", 
                      size_by = size,
                      shape_by = annoTable()$parameter
    ) + ggtitle("PCA of clustered cells") + theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
                                                  axis.text = element_text(size=12), axis.title = element_text(size=14))
    pca$guides$size$title <- "gene counts"
    
    if(length(levels(sce.qc$clusters)) <= 7 && length(unique(annoTable()$annoTable[,annoTable()$parameter]))==1) {
      pca <- pca + labs(x="PCA_1", y="PCA_2", fill="cluster") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc"))
      return(pca)
    } else if(length(levels(sce.qc$clusters)) <= 7 && length(unique(annoTable()$annoTable[,annoTable()$parameter]))>1) {
      pca <- pca + labs(x="PCA_1", y="PCA_2", colour="cluster") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc"), aesthetics = "colour")
      return(pca)
    } else {
      return(pca)
    }
  })
  
  mds_clusters <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      seuratObj <- manSeuratClustering()$seuratObj
      if(input$scaleCells_clust1 == TRUE) {
        size <- "total_features_by_counts"
      } else {
        size <- NULL
      }
    } else {
      seuratObj <- seuratClustering()$seuratObj
      if(input$scaleCells_clust == TRUE) {
        size <- "total_features_by_counts"
      } else {
        size <- NULL
      }
    }
	sce.qc <- as.SingleCellExperiment(seuratObj)
    sce.qc$clusters <- sce.qc$ident # "clusters" instead of "ident" in plot
    sce.qc <- runMDS(sce.qc, exprs_values = "logcounts")
    
    mds <- plotMDS(sce.qc, 
                   colour_by = "clusters", 
                   size_by = size,
                   shape_by = annoTable()$parameter 
    ) + ggtitle("MDS of clustered cells") + theme(plot.title = element_text(hjust = 0.5, size=18, face="plain"), 
                                                  axis.text = element_text(size=12), axis.title = element_text(size=14))
    mds$guides$size$title <- "gene counts"
    
    if(length(levels(sce.qc$clusters)) <= 7 && length(unique(annoTable()$annoTable[,annoTable()$parameter]))==1) {
      mds <- mds + labs(x="MDS_1", y="MDS_2", fill="cluster") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")) 
      return(mds)
    } else if(length(levels(sce.qc$clusters)) <= 7 && length(unique(annoTable()$annoTable[,annoTable()$parameter]))>1) {
      mds <- mds + labs(x="MDS_1", y="MDS_2", colour="cluster") + scale_fill_manual(values=c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc"), aesthetics = "colour")
      return(mds)
    } else {
      return(mds)
    }
  })
  
  ## 3d plots of clustered data
  tsne_clusters3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      seuratObj <- manSeuratClustering()$seuratObj_3D
    } else {
      seuratObj <- seuratClustering()$seuratObj_3D
    }
    
	sce.qc <- as.SingleCellExperiment(seuratObj)
    sce.qc$clusters <- sce.qc$ident # "clusters" instead of "ident" in plot
    colours <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    
    if(length(levels(sce.qc$clusters)) <= 7) {
	  plot_ly(data.frame(reducedDim(sce.qc, "TSNE", withDimnames = TRUE)), x=~tSNE_1,y=~tSNE_2, z=~tSNE_3, type="scatter3d", mode="markers", color=~sce.qc$clusters, colors=colours[1:length(levels(sce.qc$clusters))], 
              size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
        layout(sce.qc, title="3D t-SNE of clustered cells")
    } else {
	  plot_ly(data.frame(reducedDim(sce.qc, "TSNE", withDimnames = TRUE)), x=~tSNE_1,y=~tSNE_2, z=~tSNE_3, type="scatter3d", mode="markers", color=~sce.qc$clusters,
              size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
        layout(sce.qc, title="3D t-SNE of clustered cells")
    }
  })
  
  umap_clusters3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      seuratObj <- manSeuratClustering()$seuratObj_UMAP_3D
    } else {
      seuratObj <- seuratClustering()$seuratObj_UMAP_3D
    }
    
	sce.qc <- as.SingleCellExperiment(seuratObj)
    sce.qc$clusters <- sce.qc$ident # "clusters" statt "ident" im Plot
	
	UMAP_tmp <- reducedDim(sce.qc, "UMAP", withDimnames = TRUE)
	colnames(UMAP_tmp) <- c("UMAP_1", "UMAP_2", "UMAP_3")
	reducedDim(sce.qc, "UMAP", withDimnames = TRUE) <- UMAP_tmp
	rm(UMAP_tmp)
	
    colours <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    
    if(length(levels(sce.qc$clusters)) <= 7) {
	  plot_ly(data.frame(reducedDim(sce.qc, "UMAP", withDimnames = TRUE)), x=~UMAP_1,y=~UMAP_2, z=~UMAP_3, type="scatter3d", mode="markers", color=~sce.qc$clusters, colors=colours[1:length(levels(sce.qc$clusters))], 
              size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
        layout(sce.qc, title="3D UMAP of clustered cells")
    } else {
	  plot_ly(data.frame(reducedDim(sce.qc, "UMAP", withDimnames = TRUE)), x=~UMAP_1,y=~UMAP_2, z=~UMAP_3, type="scatter3d", mode="markers", color=~sce.qc$clusters, 
              size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
        layout(sce.qc, title="3D UMAP of clustered cells")
    }
  })
  
  pca_clusters3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      seuratObj <- manSeuratClustering()$seuratObj
    } else {
      seuratObj <- seuratClustering()$seuratObj
    }
	sce.qc <- as.SingleCellExperiment(seuratObj)
    sce.qc$clusters <- sce.qc$ident # "clusters" instead of "ident" in plot
	PCA_tmp <- reducedDim(sce.qc, "PCA", withDimnames = TRUE)
	colnames(PCA_tmp)[1:3] <- c("PC_1", "PC_2", "PC_3")
	reducedDim(sce.qc, "PCA", withDimnames = TRUE) <- PCA_tmp
	rm(PCA_tmp)
	
    colours <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    
    if(length(levels(sce.qc$clusters)) <= 7) {
	  plot_ly(data.frame(reducedDim(sce.qc, "PCA", withDimnames = TRUE)), x=~PC_1,y=~PC_2, z=~PC_3, type="scatter3d", mode="markers", color=~sce.qc$clusters, colors=colours[1:length(levels(sce.qc$clusters))], 
              size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
        layout(sce.qc, title="3D PCA of clustered cells")
    } else {
	  plot_ly(data.frame(reducedDim(sce.qc, "PCA", withDimnames = TRUE)), x=~PC_1,y=~PC_2, z=~PC_3, type="scatter3d", mode="markers", color=~sce.qc$clusters, 
              size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
        layout(sce.qc, title="3D PCA of clustered cells")
    }
  })
  
  mds_clusters3d <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Manual") {
      seuratObj <- manSeuratClustering()$seuratObj
    } else {
      seuratObj <- seuratClustering()$seuratObj
    }

	sce.qc <- as.SingleCellExperiment(seuratObj)
    sce.qc$clusters <- sce.qc$ident # "clusters" instead of "ident" in plot
    sce.qc <- runMDS(sce.qc, ncomponents = 3, exprs_values = "logcounts")
	
    MDS_tmp <- reducedDim(sce.qc, "MDS", withDimnames = TRUE)
	colnames(MDS_tmp) <- c("MDS_1", "MDS_2", "MDS_3")
	reducedDim(sce.qc, "MDS", withDimnames = TRUE) <- MDS_tmp
	rm(MDS_tmp)
	
    colours <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    
    if(length(levels(sce.qc$clusters)) <= 7) {
	  plot_ly(data.frame(reducedDim(sce.qc, "MDS", withDimnames = TRUE)), x=~MDS_1,y=~MDS_2, z=~MDS_3, type="scatter3d", mode="markers", color=~sce.qc$clusters, colors=colours[1:length(levels(sce.qc$clusters))], 
              size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
        layout(sce.qc, title="3D MDS of clustered cells")
    } else {
	  plot_ly(data.frame(reducedDim(sce.qc, "MDS", withDimnames = TRUE)), x=~MDS_1,y=~MDS_2, z=~MDS_3, type="scatter3d", mode="markers", color=~sce.qc$clusters, 
              size=~sce.qc$total_features_by_counts, sizes=c(20,100), text = ~paste('Cell_ID:', rownames(colData(sce.qc)), '<br>Gene count:', sce.qc$total_features_by_counts, '<br>UMI count:', sce.qc$total_counts)) %>%
        layout(sce.qc, title="3D MDS of clustered cells")
    }
  })
  
  ## heatmap of top5 marker genes of all clusters (automatic and manual analyses)
  top5_heatmap <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      markers <- seuratClustering()$markers
      seuratObj <- seuratClustering()$seuratObj
      
      # extract top5 genes per cluster
      top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
    } else if(input$automan=="Manual") {
      if(is.null(seuratMarkers()$markers) || is.null(manSeuratClustering()$seuratObj)) {
        return()
      }
      seuratObj <- manSeuratClustering()$seuratObj
      markers <- seuratMarkers()$markers
      top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
    }
    
	DoHeatmap(
      object = seuratObj,
      features = top5$gene)
  })
  
  ## heatmap of top20 marker genes of specific cluster (manual analysis)
  top20_heatmap <- reactive({
    if(is.null(seuratMarkers()$markers1) || is.null(manSeuratClustering()$seuratObj)) {
      return()
    }
    seuratObj <- manSeuratClustering()$seuratObj
    markers <- seuratMarkers()$markers1
    markers$gene <- rownames(markers)
    
    # extract top 20 genes 
    top20 <- markers %>% top_n(20, avg_logFC) 
	
    DoHeatmap(
      object = seuratObj,
      features = top20$gene)	
  })

  ## summary of all adjustable parameters
  summary <- reactive({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      sce <- build_and_filter()$sce
      seuratObj <- seuratClustering()$seuratObj
      summary_values <- c(length(sce$valid_cells), length(which(sce$valid_cells)), length(rowData(sce)$valid_cells),
                          length(which(rowData(sce)$valid_cells)), paste(">", quantile(sce$total_counts, 0.05, na.rm=TRUE)),
                          paste(">", quantile(sce$total_features_by_counts, 0.05, na.rm=TRUE)), findSignificantPCs()$elbowp, 0.8,
						  nlevels(Idents(object = seuratObj)), 0.01, 0.25)
    } else if(input$automan=="Manual") {
      sce <- geneFilter()$sce
      seuratObj <- manSeuratClustering()$seuratObj
      summary_values <- c(length(sce$valid_cells), length(which(sce$valid_cells)), length(rowData(sce)$valid_cells),
                          length(which(rowData(sce)$valid_cells)), paste(">", input$umiCutoff), paste(">", input$geneCutoff),
						  input$numPCs, input$numResolution, nlevels(Idents(object = seuratObj)), input$numPval, input$numLogfc)
    }
    summaryTable <- as.data.frame(cbind((as.character(summary_values))))
    rownames(summaryTable) <- c("Cells total", "Cells after filtering", "Genes total", "Genes after filtering", "Cutoff value UMIs",
                                "Cutoff value genes", "Significant PCs", "Resolution Clustering", "Estimated Clusters", 
                                "p-value Markergenes", "log2 fold change Markergenes")
    colnames(summaryTable) <- "Count"
    return(summaryTable)
  })
  

  ### outputs
  
  ## summary table
  output$summary <- renderTable({
    summaryTable <- summary()
    return(summaryTable)
  }, rowname=TRUE, colname=FALSE)
  
  ## histograms
  output$plot1 <- renderPlot({
    if(is.null(histo_tc())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(histo_tc())
    })
  })
  
  output$plot2 <- renderPlot({
    if(is.null(histo_tf())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(histo_tf())
    })
  })
  
  ## scatter plot genes
  output$plot2_2 <- renderPlot({
    if(is.null(scatter_genes())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(scatter_genes())
    })
  })
  
  ## raw data
  output$plot3 <- renderPlot({
    if(is.null(tsne_raw())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(tsne_raw())
    })
  })
  
  output$plot3_1 <- renderPlot({
    if(is.null(umap_raw())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(umap_raw())
    })
  })
  
  output$plot3_2 <- renderPlot({
    if(is.null(pca_raw())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(pca_raw())
    })
  })
  
  output$plot3_3 <- renderPlot({
    if(is.null(mds_raw())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(mds_raw())
    })
  })
  
  output$plot3_4 <- renderPlotly({
    if(is.null(tsne_raw3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      tsne_raw3d()
    })
  })
  
  output$plot3_5 <- renderPlotly({
    if(is.null(umap_raw3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      umap_raw3d()
    })
  })
  
  output$plot3_6 <- renderPlotly({
    if(is.null(pca_raw3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      pca_raw3d()
    })
  })
  
  output$plot3_7 <- renderPlotly({
    if(is.null(mds_raw3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      mds_raw3d()
    })
  })
  
  ## filtered data
  output$plot4 <- renderPlot({
    if(is.null(tsne_validcells())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(tsne_validcells())
    })
  })
  
  output$plot4_1 <- renderPlot({
    if(is.null(umap_validcells())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(umap_validcells())
    })
  })
  
  output$plot4_2 <- renderPlot({
    if(is.null(pca_validcells())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(pca_validcells())
    })
  })
  
  output$plot4_3 <- renderPlot({
    if(is.null(mds_validcells())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(mds_validcells())
    })
  })
  
  output$plot4_4 <- renderPlotly({
    if(is.null(tsne_validcells3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      tsne_validcells3d()
    })
  })
  
  output$plot4_5 <- renderPlotly({
    if(is.null(umap_validcells3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      umap_validcells3d()
    })
  })
  
  output$plot4_6 <- renderPlotly({
    if(is.null(pca_validcells3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      pca_validcells3d()
    })
  })
  
  output$plot4_7 <- renderPlotly({
    if(is.null(mds_validcells3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      mds_validcells3d()
    })
  })
  
  ## normalized data
  output$plot5 <- renderPlot({
    if(is.null(tsne_norm())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(tsne_norm())
    })
  })
  
  output$plot5_1 <- renderPlot({
    if(is.null(umap_norm())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(umap_norm())
    })
  })
  
  output$plot5_2 <- renderPlot({
    if(is.null(pca_norm())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(pca_norm())
    })
  })
  
  output$plot5_3 <- renderPlot({
    if(is.null(mds_norm())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(mds_norm())
    })
  })
  
  output$plot5_4 <- renderPlotly({
    if(is.null(tsne_norm3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      tsne_norm3d()
    })
  })
  
  output$plot5_5 <- renderPlotly({
    if(is.null(umap_norm3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      umap_norm3d()
    })
  })
  
  output$plot5_6 <- renderPlotly({
    if(is.null(pca_norm3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      pca_norm3d()
    })
  })
  
  output$plot5_7 <- renderPlotly({
    if(is.null(mds_norm3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      mds_norm3d()
    })
  })
  
  ## plots of normalization
  output$plot6 <- renderPlot({
    if(is.null(seuratAnalysis())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
	  return(VlnPlot_2(object = seuratAnalysis(), features.plot = c("total_features_by_counts","total_counts"), xlab = "", nCol = 2, colours = "#3c8dbc", qc = TRUE))
    })
  })
  
  output$plot7 <- renderPlot({
    if(is.null(seuratAnalysis())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(elbow())    
    })
  })
  
  output$plot8 <- renderPlot({
    if(is.null(seuratAnalysis())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      return(JackStrawPlot(object=seuratAnalysis(), dims = 1:15))
    })
  })
  
  output$plot8_2 <- renderPlot({
    if(is.null(seuratAnalysis())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
	  return(VizDimLoadings(object = seuratAnalysis(), dims=1:12, nfeatures = 10, ncol=4, reduction = "pca"))
    })
    
  })
  
  output$plot8_3 <- renderPlot({
    if(is.null(seuratAnalysis())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
	  return(DimHeatmap(object = seuratAnalysis(), dims = 1:15, nfeatures=10)) #col.use = c("#2e6c8f", "#3c8dbc", "#bc3c4d", "#8f2e3b")
    })
  })

  ## clustered data
  output$plot9 <- renderPlot({
    if(is.null(tsne_clusters())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(tsne_clusters())
    })
  })
  
  output$plot9_1 <- renderPlot({
    if(is.null(umap_clusters())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(umap_clusters())
    })
  })
  
  output$plot9_2 <- renderPlot({
    if(is.null(pca_clusters())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(pca_clusters())
    })
  })
  
  output$plot9_3 <- renderPlot({
    if(is.null(mds_clusters())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(mds_clusters())
    })
  })
  
  output$plot9_4 <- renderPlotly({
    if(is.null(tsne_clusters3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      tsne_clusters3d()
    })
  })
  
  output$plot9_5 <- renderPlotly({
    if(is.null(umap_clusters3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      umap_clusters3d()
    })
  })
  
  output$plot9_6 <- renderPlotly({
    if(is.null(pca_clusters3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      pca_clusters3d()
    })
  })
  
  output$plot9_7 <- renderPlotly({
    if(is.null(mds_clusters3d())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      mds_clusters3d()
    })
  })
  
  ## heatmap top5 marker genes of all clusters
  output$plot10 <- renderPlot({
    if(is.null(top5_heatmap())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(top5_heatmap())
    })
  })
  
  ## dropdown menu for marker genes of all clusters
  output$ui.markers <- renderUI({
    markers <- seuratClustering()$markers
    selectizeInput("markergenes", "Markergenes (max. 3):", 
                multiple = TRUE,
                options = list(maxItems = 3),
                choices= c((as.character(markers[order(markers$avg_logFC, decreasing = TRUE),]$gene))),
                selected = c(markers[order(markers$avg_logFC, decreasing = TRUE),]$gene[1:3]) 
                )
  })
  
  ## feature plot
  output$plot11 <- renderPlot({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      seuratObj <- seuratClustering()$seuratObj
      seuratObj_UMAP <- seuratClustering()$seuratObj_UMAP
      markers <- input$markergenes
    } else {
      seuratObj <- manSeuratClustering()$seuratObj
      seuratObj_UMAP <- manSeuratClustering()$seuratObj_UMAP
      markers <- input$markergenes1
    }
    
    if(input$limits=="Absolute") {
	  limits <- c(min(GetAssayData(object = seuratObj)), max(GetAssayData(object = seuratObj)))
    } else {
      limits <- c()
    }
    
    if(input$display=="t-SNE") {
      reduction <- "tsne"
      name <- "t-SNE"
    } else {
      reduction <- "umap"
      name <- "UMAP"
      seuratObj <- seuratObj_UMAP
    }
    
    if(is.null(markers)) {
      return()
    } 
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
	  fp <- FeaturePlot(
        object = seuratObj,
		features = markers,
		cols=c("lightgrey", "#3c4cbd"),
		ncol=length(markers),
		label.size = 4,
		pt.size = 2.5,
		reduction = reduction,
		combine = FALSE
      )
      fps <- lapply(fp, function(x) {x + labs(x=paste(name, "1", sep="_"), y=paste(name, "2", sep="_")) + 
          scale_colour_gradientn(limits = limits, colours=c("#2e6c8f", "#3c8dbc", "#bc3c4d", "#8f2e3b"), name="gene exp") +
          theme(plot.title=element_text(face="plain", size=18), 
                axis.text=element_text(size=12), 
                axis.title=element_text(face="plain", size=14),
                legend.title=element_text(face="plain", size=12))
      })
      cowplot::plot_grid(plotlist = fps, ncol = length(markers))
    })
  })
  
  ## violin plot of marker genes
  output$plot12 <- renderPlot({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      seuratObj <- seuratClustering()$seuratObj
      markers <- input$markergenes
    } else {
      seuratObj <- manSeuratClustering()$seuratObj
      markers <- input$markergenes1
    }
    
    if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      VlnPlot_2(object = seuratObj, features.plot = c(markers), xlab = "Cluster", nCol = length(markers), colours = colour, qc = FALSE)
    })
  })
  
  ## ridge plot
  output$plot12_2 <- renderPlot({
    if(is.null(input$automan)) {
      return()
    } else if(input$automan=="Automatic") {
      seuratObj <- seuratClustering()$seuratObj
      markers <- input$markergenes
    } else {
      seuratObj <- manSeuratClustering()$seuratObj
      markers <- input$markergenes1
    }
    
	if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      RidgePlot_2(object = seuratObj, features.plot = c(markers), nCol = length(markers), colours = colour)
    })    
  })
  
  ## dot plot
  output$plot12_3 <- renderPlot({
    if(is.null(dotplot1())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(dotplot1())
    })
  })
  
  ## heatmap with selected genes
  output$plot10_2 <- renderPlot({
    if(is.null(heatMark1())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(heatMark1())
    })
  })
  
  ## feature plot (co-expression)
  output$plot11_2 <- renderPlot({
    if(is.null(coEx1())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      coEx1()
    })
  })
  
  ## heatmap top20 marker genes of a specific cluster
  output$plot13 <- renderPlot({
    if(is.null(top20_heatmap())) {
      return()
    }
    withProgress(message="Clustering Cells.", detail="This might take some time.", value = 1, {
      plot(top20_heatmap())
    })
  })
  
  ## feature plot
  output$plot14 <- renderPlot({
    seuratObj <- manSeuratClustering()$seuratObj
    seuratObj_UMAP <- manSeuratClustering()$seuratObj_UMAP
    markers <- input$markergenes2
    
    if(input$limits2=="Absolute") {
	  limits <- c(min(GetAssayData(object = seuratObj)), max(GetAssayData(object = seuratObj)))
    } else {
      limits <- c()
    }
    
    if(input$display2=="t-SNE") {
      reduction <- "tsne"
      name <- "t-SNE"
    } else {
      reduction <- "umap"
      name <- "UMAP"
      seuratObj <- seuratObj_UMAP
    }
    
    if(is.null(markers)) {
      return()
    } else {
      withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
		fp <- FeaturePlot(
          object = seuratObj,
		  features = markers,
		  cols=c("lightgrey", "#3c4cbd"),
		  ncol=length(markers),
		  label.size = 4,
		  pt.size = 2.5,
		  reduction = reduction,
		  combine = FALSE
        )
        fps <- lapply(fp, function(x) {x + labs(x=paste(name, "1", sep="_"), y=paste(name, "2", sep="_")) + 
            scale_colour_gradientn(limits = limits, colours=c("#2e6c8f", "#3c8dbc", "#bc3c4d", "#8f2e3b"), name="gene exp") +
            theme(plot.title=element_text(face="plain", size=18), 
                  axis.text=element_text(size=12), 
                  axis.title=element_text(face="plain", size=14),
                  legend.title=element_text(face="plain", size=12))
        })
        cowplot::plot_grid(plotlist = fps, ncol = length(markers))
      })
    }
  })
  
  ## violin plot of marker genes
  output$plot15 <- renderPlot({
    seuratObj <- manSeuratClustering()$seuratObj
    markers <- input$markergenes2
	if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      VlnPlot_2(object=seuratObj, features.plot = c(markers), xlab = "Cluster", nCol = length(markers), colours = colour, qc = FALSE)
    })
  })
  
  ## ridge plot
  output$plot15_2 <- renderPlot({
    seuratObj <- manSeuratClustering()$seuratObj
    markers <- input$markergenes2
    
	if(length(levels(seuratObj)) <= 7) {
      colour <- c("#bc6b3c", "#3c8dbc", "#bc3c8d", "#3cbc6d", "#3c4dbc", "#bc3c4d", "#6b3cbc")
    } else {
      colour <- NULL
    }
    
    if(is.null(markers)) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      RidgePlot_2(object=seuratObj, features.plot = c(markers), nCol = length(markers), colours = colour)
    })
  })
  
  ## dot plot
  output$plot15_3 <- renderPlot({
    if(is.null(dotplot2())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(dotplot2())
    })
  })
  
  ## heatmap with selected genes
  output$plot13_2 <- renderPlot({
    if(is.null(heatMark2())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      plot(heatMark2())
    })
  })
  
  ## feature plot (co-expression)
  output$plot14_2 <- renderPlot({
    if(is.null(coEx2())) {
      return()
    }
    withProgress(message="Creating Plot.", detail="This might take some time.", value = 1, {
      coEx2()
    })
  })
  
  
  ### contact
  
  output$imp <- renderText({
    paste0("<p>wasp@computational.bio</p>
	        <p>andreas.hoek@computational.bio.uni-giessen.de</p><br>",
			"<H2>Web:</H2>",
           "<p><a href='http://www.computational.bio/software/wasp/' target='_blank'>http://www.computational.bio/software/wasp/</a></p>
		    <p><a href='https://www.uni-giessen.de/fbz/fb08/Inst/bioinformatik' target='_blank'>https://www.uni-giessen.de/fbz/fb08/Inst/bioinformatik</a></p>")
  })

  
  ### downloads
  
  ## function for customizing the feature plot of Seurat
  FeaturePlot_2 <- function(fp, sObj, limits, name) {
    main_function <- function(fp=fp, sObj1=sObj, limits1=limits, name1=name) {
      fp <- fp + labs(x=paste(name1, "1", sep="_"), y=paste(name1, "2", sep="_")) + 
        scale_colour_gradientn(limits = limits1, colours=c("#2e6c8f", "#3c8dbc", "#bc3c4d", "#8f2e3b"), name="gene exp") +
        theme(plot.title=element_text(face="plain", size=18), 
              axis.text=element_text(size=12), 
              axis.title=element_text(face="plain", size=14),
              legend.title=element_text(face="plain", size=12))
    }
    
    fp <- lapply(X=fp, FUN=main_function)
    
    return(cowplot::plot_grid(plotlist = fp, ncol = 1))
  }
  
  ## feature plot for download (automatic analysis)
  FeaturePlotAUTO <- reactive({
    seuratObj <- seuratClustering()$seuratObj
    seuratObj_UMAP <- seuratClustering()$seuratObj_UMAP
    markers <- input$markergenes
    
    if(input$display=="t-SNE") {
      reduction <- "tsne"
      name <- "t-SNE"
    } else {
      reduction <- "umap"
      name <- "UMAP"
      seuratObj <- seuratObj_UMAP
    }
    
    if(is.null(markers)) {
      return()
    }
    fp <- FeaturePlot(object = seuratObj, features = markers,
                          ncol=1, label.size = 4, pt.size = 1.5, reduction = reduction, combine = FALSE)						  
  })
  
  ## feature plot for download (manual analysis) (all clusters)
  FeaturePlotMAN1 <- reactive({
    seuratObj <- manSeuratClustering()$seuratObj
    seuratObj_UMAP <- manSeuratClustering()$seuratObj_UMAP
    markers <- input$markergenes1
    
    if(input$display=="t-SNE") {
      reduction <- "tsne"
      name <- "t-SNE"
    } else {
      reduction <- "umap"
      name <- "UMAP"
      seuratObj <- seuratObj_UMAP
    }
    
    if(is.null(markers)) {
      return()
    }
    fp <- FeaturePlot(object = seuratObj, features = markers,
                          ncol=1, label.size = 4, pt.size = 1.5, reduction = reduction, combine = FALSE)
  })
  
  ## feature plot for download (manual analysis) (specific cluster)
  FeaturePlotMAN2 <- reactive({
    seuratObj <- manSeuratClustering()$seuratObj
    seuratObj_UMAP <- manSeuratClustering()$seuratObj_UMAP
    markers <- input$markergenes2
    
    if(input$display2=="t-SNE") {
      reduction <- "tsne"
      name <- "t-SNE"
    } else {
      reduction <- "umap"
      name <- "UMAP"
      seuratObj <- seuratObj_UMAP
    }
    
    if(is.null(markers)) {
      return()
    }
						  
	fp <- FeaturePlot(object = seuratObj, features = markers,
                          ncol=1, label.size = 4, pt.size = 1.5, reduction = reduction, combine = FALSE)
  })

  ## download all plots
  output$report <- downloadHandler(
    filename = "plots.pdf",
    content = function(files) {
      pdf(files)
      if(is.null(input$automan)) {
        return()
      } else if(input$automan=="Automatic") {
        plot(histo_tc())
        plot(histo_tf())
        plot(scatter_genes())
        plot(tsne_raw())
        plot(umap_raw())
        plot(pca_raw())
        plot(mds_raw())
        plot(tsne_validcells())
        plot(umap_validcells())
        plot(pca_validcells())
        plot(mds_validcells())
        plot(tsne_norm())
        plot(umap_norm())
        plot(pca_norm())
        plot(mds_norm())
        plot(violin())
        plot(elbow())
        print(JackStrawPlot(object=seuratAnalysis(), dims = 1:15))
		print(VizDimLoadings(object = seuratAnalysis(), dims=1:12, nfeatures = 10, ncol=4, reduction = "pca"))
		print(DimHeatmap(object = seuratAnalysis(), dims = 1:15, nfeatures=10))
        plot(tsne_clusters())
        plot(umap_clusters())
        plot(pca_clusters())
        plot(mds_clusters())
        plot(top5_heatmap())
        
        if(!is.null(input$markergenes)) {
          if(input$limits=="Absolute") {
            if(input$display=="t-SNE") {
              print(FeaturePlot_2(fp=FeaturePlotAUTO(), sObj=seuratClustering()$seuratObj, limits=c(min(GetAssayData(seuratClustering()$seuratObj)), max(GetAssayData(seuratClustering()$seuratObj))), name="t-SNE"))
            } else {
			  print(FeaturePlot_2(fp=FeaturePlotAUTO(), sObj=seuratClustering()$seuratObj_UMAP, limits=c(min(GetAssayData(seuratClustering()$seuratObj_UMAP)), max(GetAssayData(seuratClustering()$seuratObj_UMAP))), name="UMAP"))
            }
          } else {
            if(input$display=="t-SNE") {
              print(FeaturePlot_2(fp=FeaturePlotAUTO(), sObj=seuratClustering()$seuratObj, limits=c(), name="t-SNE"))
            } else {
              print(FeaturePlot_2(fp=FeaturePlotAUTO(), sObj=seuratClustering()$seuratObj_UMAP, limits=c(), name="UMAP"))
            }
          }
          plot(violin1())
          plot(ridge1())
          plot(dotplot1())
          plot(heatMark1())
          print(coEx1())
        }
      } else {
        plot(histo_tc())
        plot(histo_tf())
        plot(scatter_genes())
		print("scatter")
        plot(tsne_raw())
        plot(umap_raw())
        plot(pca_raw())
        plot(mds_raw())
		print("raw")
        plot(tsne_validcells())
        plot(umap_validcells())
        plot(pca_validcells())
        plot(mds_validcells())
		print("valid_cells")
        plot(tsne_norm())
        plot(umap_norm())
        plot(pca_norm())
        plot(mds_norm())
		print("normalized_data")
        plot(violin())
		print("violin")
        plot(elbow())
        print(JackStrawPlot(object=seuratAnalysis(), dims = 1:15)) 
		print("normalization_pc_selection")
		print(VizDimLoadings(object = seuratAnalysis(), dims=1:12, nfeatures = 10, ncol=4, reduction = "pca"))
		print(DimHeatmap(object = seuratAnalysis(), dims = 1:15, nfeatures=10))
		print("heatmap")
        plot(tsne_clusters())
        plot(umap_clusters())
        plot(pca_clusters())
        plot(mds_clusters())
		print("clusters")
        if(input$allMarkers && input$find!="None") {
          plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
          text(2.3, 2.6, "Markergenes of all Clusters", pos=1, cex=2.5)
          plot(top5_heatmap())
          
          if(!is.null(input$markergenes1)) {
            if(input$limits=="Absolute") {
              if(input$display=="t-SNE") {
				print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj)), max(GetAssayData(manSeuratClustering()$seuratObj))), name="t-SNE"))
              } else {
                print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj_UMAP)), max(GetAssayData(manSeuratClustering()$seuratObj_UMAP))), name="UMAP"))
              }
            } else {
              if(input$display=="t-SNE") {
                print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj, limits=c(), name="t-SNE"))
              } else {
                print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(), name="UMAP"))
              }
            }
            
            plot(violin1())
            plot(ridge1())
            plot(dotplot1())
            plot(heatMark1())
            print(coEx1())
            
          }
          plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
          text(2.3, 2.6, paste0("Markergenes of Cluster", input$find, sep=""), pos=1, cex=2.5)
          plot(top20_heatmap())
          
          if(!is.null(input$markergenes2)) {
            if(input$limits2=="Absolute") {
              if(input$display2=="t-SNE") {
                print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj)), max(GetAssayData(manSeuratClustering()$seuratObj))), name="t-SNE"))
              } else {
				print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj_UMAP)), max(GetAssayData(manSeuratClustering()$seuratObj_UMAP))), name="UMAP"))
              }
            } else {
              if(input$display2=="t-SNE") {
                print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj, limits=c(), name="t-SNE"))
              } else {
                print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(), name="UMAP"))
              }
            }
            
            plot(violin2())
            plot(ridge2())
            plot(dotplot2())
            plot(heatMark2())
            print(coEx2())
          }
        } else if(input$allMarkers && input$find=="None") {
          plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
          text(2.3, 2.6, "Markergenes of all Clusters", pos=1, cex=2.5)
          plot(top5_heatmap())
          
          if(!is.null(input$markergenes1)) {
            if(input$limits=="Absolute") {
              if(input$display=="t-SNE") {
                print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj, limits=c(min(manSeuratClustering()$seuratObj@data), max(manSeuratClustering()$seuratObj@data)), name="t-SNE"))
              } else {
                print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(min(manSeuratClustering()$seuratObj_UMAP@data), max(manSeuratClustering()$seuratObj_UMAP@data)), name="UMAP"))
              }
            } else {
              if(input$display=="t-SNE") {
                print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj, limits=c(), name="t-SNE"))
              } else {
                print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(), name="UMAP"))
              }
            }
            
            plot(violin1())
            plot(ridge1())
            plot(dotplot1())
            plot(heatMark1())
            print(coEx1())
          }
        } else if(input$allMarkers==FALSE && input$find!="None") {
          plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
          text(2.3, 2.6, paste0("Markergenes of Cluster", input$find, sep=""), pos=1, cex=2.5)
          plot(top20_heatmap())
          
          if(!is.null(input$markergenes2)) {
            if(input$limits2=="Absolute") {
              if(input$display2=="t-SNE") {
				print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj)), max(GetAssayData(manSeuratClustering()$seuratObj))), name="t-SNE"))
              } else {
				print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj_UMAP)), max(GetAssayData(manSeuratClustering()$seuratObj_UMAP))), name="UMAP"))
              }
            } else {
              if(input$display2=="t-SNE") {
                print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj, limits=c(), name="t-SNE"))
              } else {
                print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(), name="UMAP"))
              }
            }
            
            plot(violin2())
            plot(ridge2())
            plot(dotplot2())
            plot(heatMark2())
            print(coEx2())
          }
        } else {
          
        }
      }
      dev.off()
    }
  )
  
  ## download plots of filtering (manual analysis)
  output$filter_man <- downloadHandler(
    filename = "filtering_plots.pdf",
    content = function(files) {
      pdf(files)
      plot(histo_tc())
      plot(histo_tf())
      plot(scatter_genes())
      dev.off()
    }
  )
  
  ## dowload plots of normalization (manual analysis)
  output$norm_man <- downloadHandler(
    filename = "normalisation_plots.pdf",
    content = function(files) {
      pdf(files)
      print(elbow())
      print(JackStrawPlot(object=seuratAnalysis(), dims = 1:15))
	  print(VizDimLoadings(object = seuratAnalysis(), dims=1:12, nfeatures = 10, ncol=4, reduction = "pca"))
	  print(DimHeatmap(object = seuratAnalysis(), dims = 1:15, nfeatures=10))
      dev.off()
    }
  )
  
  ## download plots of clustering (manual analysis)
  output$cluster_man <- downloadHandler(
    filename = "clustering_plot.pdf",
    content = function(files) {
      pdf(files)
      plot(tsne_clusters())
      plot(umap_clusters())
      plot(pca_clusters())
      plot(mds_clusters())
      dev.off()
    }
  )
  
  ## download plots of marker gene calculation (manual analysis)
  output$marker_man <- downloadHandler(
    filename = "markergenes_plots.pdf",
    content = function(files) {
      pdf(files)
      if(input$allMarkers && input$find!="None") {
        plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
        text(2.3, 2.6, "Markergenes of all Clusters", pos=1, cex=2.5)
        plot(top5_heatmap())
        
        if(!is.null(input$markergenes1)) {
          if(input$limits=="Absolute") {
            if(input$display=="t-SNE") {
			  print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj)), max(GetAssayData(manSeuratClustering()$seuratObj))), name="t-SNE"))
            } else {
			  print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj_UMAP)), max(GetAssayData(manSeuratClustering()$seuratObj_UMAP))), name="UMAP"))
            }
          } else {
            if(input$display=="t-SNE") {
			  print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj, limits=c(), name="t-SNE"))
            } else {
              print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(), name="UMAP"))
            }
          }
          
          plot(violin1())
          plot(ridge1())
          plot(dotplot1())
          plot(heatMark1())
          print(coEx1())
        }
        plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
        text(2.3, 2.6, paste0("Markergenes of Cluster", input$find, sep=""), pos=1, cex=2.5)
        plot(top20_heatmap())
        
        if(!is.null(input$markergenes2)) {
          if(input$limits2=="Absolute") {
            if(input$display2=="t-SNE") {
			  print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj)), max(GetAssayData(manSeuratClustering()$seuratObj))), name="t-SNE"))
            } else {
              print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj_UMAP)), max(GetAssayData(manSeuratClustering()$seuratObj_UMAP))), name="UMAP"))
            }
          } else {
            if(input$display2=="t-SNE") {
              print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj, limits=c(), name="t-SNE"))
            } else {
              print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(), name="UMAP"))
            }
          }
          
          plot(violin2())
          plot(ridge2())
          plot(dotplot2())
          plot(heatMark2())
          print(coEx2())
        }
      } else if(input$allMarkers && input$find=="None") {
        plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
        text(2.3, 2.6, "Markergenes of all Clusters", pos=1, cex=2.5)
        plot(top5_heatmap())
       
        if(!is.null(input$markergenes1)) {
          if(input$limits=="Absolute") {
            if(input$display=="t-SNE") {
			  print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj)), max(GetAssayDatamanSeuratClustering()$seuratObj)), name="t-SNE"))
            } else {
			  print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj_UMAP)), max(GetAssayData(manSeuratClustering()$seuratObj_UMAP))), name="UMAP"))
            }
          } else {
            if(input$display=="t-SNE") {
              print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj, limits=c(), name="t-SNE"))
            } else {
              print(FeaturePlot_2(fp=FeaturePlotMAN1(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(), name="UMAP"))
            }
          }
          
          plot(violin1())
          plot(ridge1())
          plot(dotplot1())
          plot(heatMark1())
          print(coEx1())
        }
      } else if(input$allMarkers==FALSE && input$find!="None"){
        plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
        text(2.3, 2.6, paste0("Markergenes of Cluster", input$find, sep=""), pos=1, cex=2.5)
        plot(top20_heatmap())
        
        if(!is.null(input$markergenes2)) {
          if(input$limits2=="Absolute") {
            if(input$display2=="t-SNE") {
			  print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj)), max(GetAssayData(manSeuratClustering()$seuratObj))), name="t-SNE"))
            } else {
			  print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(min(GetAssayData(manSeuratClustering()$seuratObj_UMAP)), max(GetAssayData(manSeuratClustering()$seuratObj_UMAP))), name="UMAP"))
            }
          } else {
            if(input$display2=="t-SNE") {
              print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj, limits=c(), name="t-SNE"))
            } else {
              print(FeaturePlot_2(fp=FeaturePlotMAN2(), sObj=manSeuratClustering()$seuratObj_UMAP, limits=c(), name="UMAP"))
            }
          }
          
          plot(violin2())
          plot(ridge2())
          plot(dotplot2())
          plot(heatMark2())
          print(coEx2())
        }
      } else {
        return()
      }
      dev.off()
    }
  )
  
  ## download normalization table
  output$normtable <- downloadHandler(
    filename = "normalisation_table.csv",
    content = function(file) {
      write.table(normMatrix(), file, sep=",", row.names=TRUE, col.names=NA)
    }
  )
  
  ## download normalization table (manual analysis)
  output$normtable_man <- downloadHandler(
    filename = "normalisation_table.csv",
    content = function(file) {
      write.table(normMatrix(), file, sep=",", row.names=TRUE, col.names=NA)
    }
  )
  
  ## download table of marker genes (automatic analysis)
  output$markerlist <- downloadHandler(
    filename = "markergenes.csv",
    content = function(file) {
      write.table(seuratClustering()$markers, file, sep=",", row.names=TRUE, col.names=NA)
    }
  )
  
  ## download table of marker genes (all clusters) (manual analysis)
  output$downloadAllMarkers <- downloadHandler(
    filename = "markergenes_allClusters.csv",
    content = function(file) {
      write.table(seuratMarkers()$markers, file, sep=",", row.names=TRUE, col.names=NA)
    }
  )

  output$downloadAllMarkers1 <- downloadHandler(
    filename = "markergenes_allClusters.csv",
    content = function(file) {
      write.table(seuratMarkers()$markers, file, sep=",", row.names=TRUE, col.names=NA)
    }
  )
  
  ## download table of marker genes (specific cluster) (manual analysis)
  output$downloadSpecMarkers <- downloadHandler(
    filename = paste("margergenes_cluster", input$find, ".csv", sep=""),
    content = function(file) {
      markers <- seuratMarkers()$markers1
      markers$gene <- rownames(markers)
      write.table(markers, file, sep=",", row.names=TRUE, col.names=NA)
    }
  )
  
  output$downloadSpecMarkers1 <- downloadHandler(
    filename = paste("margergenes_cluster", input$find, ".csv", sep=""),
    content = function(file) {
      markers <- seuratMarkers()$markers1
      markers$gene <- rownames(markers)
      write.table(markers, file, sep=",", row.names=TRUE, col.names=NA)
    }
  )
  
  ## download summary
  output$logfile <- downloadHandler(
    filename = "summary.csv",
    content = function(file) {
      write.table(summary(), file, sep=",", row.names=TRUE, col.names=FALSE)
    }
  )
  
  ## stop R script if website is closed
  session$onSessionEnded(stopApp) 
}

shinyApp(ui = ui, server = server)