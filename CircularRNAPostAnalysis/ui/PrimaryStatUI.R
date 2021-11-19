
library(highcharter)

tabItem(tabName = "PrimaryStat",
        #busyIndicator("Calculation In progress", wait = 0),
        headerPanel("Primary Statistics"),
        fluidRow(
                 column(6, 
                     box(title = "Classes distributions Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12, 
                         actionButton("PrintPlot", "Display Class plot"), 
                         highchartOutput("ClassPlot",height = "500px")
                     )
                 ),   
                 column(6,
                        box(title = "Classes distributions Table", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12,
                          DT::dataTableOutput("ClassTable")
                        )
                      )
                 
          ),
        fluidRow(
          column(6,
                 box(title = "Length distributions Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12, 
                     actionButton("PrintPlotLength", "Display Length plot"), 
                     highchartOutput("LengthPlot",height = "500px")
                 )
          ),
          column(6,
                box(title = "Length distributions Table", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12,
                DT::dataTableOutput("LengthTable")
                )
          )
        ),
        fluidRow(
          column(6,
                 box(title = "Signal distributions Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12, 
                     actionButton("PrintPlotSignal", "Display Signal plot"), 
                     highchartOutput("SignalPlot",height = "500px")
                 )
          ),
          column(6,
                 box(title = "Signal distributions Table", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12,
                     DT::dataTableOutput("SignalTable")
                 )
          )
        ),
        fluidRow(
          column(6,
                 box(title = "Supported reads distributions Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12, 
                     actionButton("PrintPlotSuppReads", "Display Supported reads plot"), 
                     highchartOutput("SuppReadsPlot",height = "500px")
                 )
          ),
          column(6,
                 box(title = "Supported reads distributions Table", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12,
                     DT::dataTableOutput("SuppReadsTable")
                 )
          )
        ),
        fluidRow(
          column(6,
                 box(title = "Supported Method distributions Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12, 
                     actionButton("PrintPlotSuppMeth", "Display Supported methods plot"), 
                     highchartOutput("SuppMethPlot",height = "500px")
                 )
          ),
          column(6,
                 box(title = "Supported Method distributions Table", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12,
                     DT::dataTableOutput("SuppMethTable")
                 )
          )
        ),
        fluidRow(
          box(title = "Tools statistics", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, width = 12, 
              actionButton("PrintPlotToolsStatistics", "Display Tools statistics plot"), 
              highchartOutput("ToolsStatisticsPlot",height = "500px")
          )
        )
        
        
)

