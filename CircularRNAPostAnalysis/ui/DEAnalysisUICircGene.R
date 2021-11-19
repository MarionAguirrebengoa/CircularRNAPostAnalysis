tabItem(tabName = "DEOnCircGene",
        #busyIndicator("Calculation In progress", wait = 0),
        conditionalPanel(condition = "output.TestIfDeIsDoneOnGeneANDCircular", 
                         helpText("DE Analysis on gene or on circulars isn't done !"),

        headerPanel("Common DE Analysis on Circular And Gene"),
        # box(title = "Table of common Up and Down regulated between circular RNA genes counts and global genes counts", width = 12,
        #     DT::dataTableOutput("CommonUpDown")
        # ),
        box(width = 4, align = "center", status = "warning",
            actionButton("CAH_Button", "Compute CAH")
        ),
        br(),
        box(title = "Dendrogramme of CAH", width = 12,
        plotOutput("Dendrogramme_CAH")
        ),
        box(title = "Inertie of CAH", width = 12,
            plotOutput("Intertie")
        ),
        fluidRow(
            box(title = "Scatter plot CAH parameters", width = 12,
                  numericInput(inputId = "NbClasses", label = "Number of Classes", value = 2, min = 0),
                  actionButton("CAH_ScatterPlot_Button", "Load Plot"),
                  uiOutput("SelectGeneId")
            )
        ),
        box(title = "Scatter Plot CAH", width = 12,
            highchartOutput("CAH_ScatterPlot", height = "500px")
        ),
        downloadButton("download_scatter_CAH_points", "Download results from CAH (points)"),
        downloadButton("download_scatter_CAH_barycentres", "Download results from CAH (barycentres)")
        
        
  )      
)