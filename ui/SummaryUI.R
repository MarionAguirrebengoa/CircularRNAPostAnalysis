tabItem(tabName = "Summary",
        #busyIndicator("Calculation In progress", wait = 0),
        headerPanel("Summary"),
        fluidRow(column(12, align = "center",
                       infoBoxOutput("nbReplicats"),
                       infoBoxOutput("nbCdn"),
                       infoBoxOutput("nbSample")
                 ) 
        ),
        fluidRow(align = "center",
                 checkboxInput("isStrandedData", "is stranded", FALSE),
                 box(width = 12, title = "Circular RNA prediction resume", status = "info",
                       tableOutput("summary")
                 )
        ),
        fluidRow(align = "center",
                 box(width = 12, title = "Processing Informations", status = "info",
                       DT::dataTableOutput("summaryProcessing")
                 )
        ),
        fluidRow(align = "center",
                 box(width = 12, title = "Merging Informations", status = "info",
                       DT::dataTableOutput("MergingInformations")
                 )
        ),
        actionButton("ExportToRDS", "Export Table to .rds file"),
                 fluidRow(
                   box(width = 12, title = "View circular RNA predictions", status = "warning",
                           column(4,
                              textInput("dataBaseGB", label = h5("DataBase"), value = "Enter a database (e.g hg38)")
                           ),
                           column(4,
                              textInput("OrganismGB", label = h5("Organism"), value = "Enter an organism (e.g human)")
                           ),
                           column(4,
                              actionButton("GenerateLinksGB", "Generate Links To UCSC")
                           )
                      )
                 ),
                 fluidRow(
                     box(width = 12, status = "success",
                           downloadButton("download_selected_rows", "Download selected rows"),
                           DT::dataTableOutput("circRNATable")
                     )
                 )
        
        
)