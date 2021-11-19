tabItem(tabName = "DEOnCirc",
            #busyIndicator("Calculation In progress", wait = 0),
               fluidRow(
                 column(12,
                          headerPanel("Differential expression analysis on circular RNA"),
                          br(),
                          br()
                        )
               ),
          fluidRow(
            column(4,
                  box(width=12,  status = "warning",
                  actionButton("reloadCountMatrix", "Load/Update Count matrix"),
                  numericInput(inputId = "JunctionCutOff", label = "Circular Junction Cutoff", value = 0, min = 0)
                  ),
                  box(title = "Fitlers on counts", width = 12, status = "warning",
                      numericInput(inputId = "BasemeanCutCirc", label = "BaseMean Cutoff on genes", value = 0, min = 0),
                      numericInput(inputId = "BasemeanCutJunctions", label = "BaseMean Cutoff on junctions", value = 0, min = 0)
                  )
          ),
          column(4,
                 box(title = "Analyses", width = 12, status = "warning",
                       actionButton("DEAnalysisOnlyCircButt", "Differential Expression Analysis (Circular on gene)"),
                       br(),
                       actionButton("KSTestButt", "Compute Kolmogorov-Smirnov Test"),
                       br(),
                       actionButton("DEAnalysisOnlyJunctionsButt", "Differential Expression Analysis (Circular Junctions)")
                 ),
                 box(title = "Cut off for Up and Down regulated genes", width = 12, status = "warning",
                     column(6,
                            numericInput("CutOffFC_Circ", "CutOff for Log2 Fold-Change", 0.5, min = 0)
                     ),
                     column(6,
                            numericInput("CutPadj_Circ", "CutOff for p.adjust", 0.05, min = 0)
                     )
                 )
          ),
          column(4,
                 fluidRow(
                   infoBoxOutput("Total_number_of_genes", width = 12)
                 ),
                 fluidRow(
                   infoBoxOutput("Number_of_genes_Raw_counts", width = 12)
                 ),
                 fluidRow(
                   infoBoxOutput("Number_of_genes_Filt_counts", width = 12)
                 )
                 
          )
          ),
          fluidRow(
            box(title = "Normalisation method", width = 12, status = "warning",
                  column(6,
                         radioButtons("NormMethod_Circ", label = "Choose your normalisation method",
                                      choices = list("normalisation factor from genes" = 1, "CPM" = 2, "DeSeq" = 3, "RPM" = 4), 
                                      selected = 3),
                         br(),
                         helpText("CPM (Count per million) is done using the number of reads mapped only on GTF features "),
                         helpText("RPM (Read per million) is done using the number of reads aligned on all the genome (including chimeric reads)")
                  ),
                  column(6,
                        numericInput("normFactCPM_Circ", label = "Per million factor : (count/libSize) x perMillionFactor", value = 1000000, , min = 0)
                  )
            )
          ),
          br(),
          fluidRow(
            column(3, 
                   box(title = "Normalisation Factors for Circular RNAs by genes", width = 12, status = "info",
                      tableOutput("normFactorCircByGenes")
                   )
                   ),
            column(6,
                   box(title = "Model parameters", width = 12, status = "warning",
                       column(4,
                              uiOutput("mainEffectCirc"),
                              uiOutput("batchEffectCirc")
                       ),
                       column(4,
                              uiOutput("Condition_1_circ"),
                              uiOutput("Condition_2_circ")
                       ),
                       column(4,
                              radioButtons("DE_Test_circ", label = "Choose your test",
                                           choices = list("Wald" = 1, "LRT" = 2), 
                                           selected = 1)
                       )
                       
                   )
                   
                   ),
            column(3,
                   box(title = "Normalisation Factors for Circular RNAs junctions", width = 12, status = "info",
                       tableOutput("normFactorCircJunctions")
                       )
                   )
            
          ),
fluidRow(
      column(6,
             textOutput("Cdn1")
      ),
      column(6,
             textOutput("Cdn2")
      )
  ),

  box(title = "Boxplot of raw counts", width = 12, status = "success",
    highchartOutput("circDensplot",height = "500px")
  ),
  box(title = "Quantile plot of raw counts", width = 12, status = "success",
    highchartOutput("QuantilePlot",height = "500px")
  ),
  br(),
  fluidRow(
        box(title = "Dispersion plot for normalised counts (on genes)", width = 12, status = "success",
              column(6,
                     uiOutput("Rep1_Circ")
                     ),
              column(6,
                     uiOutput("Rep2_Circ")
                     ),
              highchartOutput("DispersionPlot",height = "500px"),
              br(),
              column(6,
                     textOutput("Rsquared_Cdn1_Circ")
              ),
              column(6,
                     textOutput("Rsquared_Cdn2_Circ")
              )
        ),
        br(),
        box(title = "Normalised counts filter by baseMean", width = 12, status = "success",
              highchartOutput("circDensplotNormalised",height = "500px")
        ),
        br(),
        box(title = "Differential expression analysis on Circular group by genes", width = 12, status = "success",
              DT::dataTableOutput("DEResultsTableOnlyOnCirc")
        ),
        box(title = "Table for Up and Down regulated genes", width = 12, status = "success",
            column(6,
                  highchartOutput("PlotUpDownForCirc")
            ),
            column(6,
                   DT::dataTableOutput("TableUpDownForCirc") 
                   )
        ),
      
        box(title = "Kolmogorov-Smirnov test", width = 12, status = "success",
              DT::dataTableOutput("KS_test_Table")
        ),
        box(title = "Differential expression analysis on Circular Junctions", width = 12,  status = "success",
              DT::dataTableOutput("DEResultsTableOnlyOnJunctionsCirc")
        ),
        box(title = "Table for Up and Down regulated junctions", width = 12,  status = "success",
            column(6,
                  highchartOutput("PlotUpDownForJunctionsCirc")
            ),
            column(6,
                   DT::dataTableOutput("TableUpDownForJunctionsCirc") 
            )
        ),
        
       
        box(title = "Normalised counts for Circular RNAs by genes", width = 12,  status = "success",
              DT::dataTableOutput("CountsNormalisedForCircularOnGenes")
        ),
        
        box(title = "Raw counts for circular RNA junctions building genes", width = 12,  status = "success",
              DT::dataTableOutput("RawCountsForJunctionsThatBuildGenes")
        ),
      
        box(title = "Normalised counts for Circular RNAs junctions", width = 12, status = "success",
              DT::dataTableOutput("CountsNormalisedForCircularJunctions")
        )
        )      
  
)