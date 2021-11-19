tabItem(tabName = "DEOnGene",
        #busyIndicator("Calculation In progress", wait = 0),
        headerPanel("Differential expression analysis on Gene"),
        br(),
        fluidRow(
        column(4,
               box(title = "Basemean cut off", width = 12, status = "warning",
                  numericInput(inputId = "BasemeanCutGene", label = "BaseMean Cutoff", value = 0, min = 0)
               )
               ),
        column(4,
               box(title = "Normalisation methods", width = 12, status = "warning",
                     radioButtons("NormMethod", label = "Choose your normalisation method",
                                  choices = list("House Keeping Genes" = 1, "CPM" = 2, "DeSeq" = 3, "RPM" = 4), 
                                  selected = 3)
               ),
               helpText("CPM (Count per million) is done using the number of reads mapped only on GTF features "),
               helpText("RPM (Read per million) is done using the number of reads aligned on all the genome (including chimeric reads)")
        ),
        column(4,
               box(title = "Normalisation parameters", width = 12, status = "warning",
                    uiOutput("HouseKeepingGenes"),
                    numericInput("normFactCPM", label = "Per million factor : (count/libSize) x perMillionFactor", value = 1000000, min = 0)
               )
        )
        ),
        br(),
        fluidRow
        (
          column(3,
              box(title = "Cross reference paramters", width = 12, status = "warning",
                    column(6,
                           selectInput("Idtype", "Source Ids", choices = c("Ensembl"="ensembl_gene_id", "RefSeq" = "entrezgene", "None"="None"), selected =  "None")
                    ),
                    column(6,
                           uiOutput("Organism")
                    )
              )
          ),
          column(3,
              box(title = "Cut off for Up and Down regulated genes", width = 12, status = "warning",
                column(6,
                       numericInput("CutOffFC_Gene", "CutOff for Log2 Fold-Change", 0.5, min = 0)
                ),
                column(6,
                       numericInput("CutPadj_Gene", "CutOff for p.adjust", 0.05, min = 0)
                )
              )
          ),
          column(6,
                 box(title = "Model parameters", width = 12, status = "warning",
                           column(4,
                                  uiOutput("mainEffect"),
                                  uiOutput("batchEffect")
                           ),
                           column(4,
                                  uiOutput("Condition_1_gene"),
                                  uiOutput("Condition_2_gene")
                           ),
                           column(4,
                                  radioButtons("DE_Test_gene", label = "Choose your test",
                                              choices = list("Wald" = 1, "LRT" = 2), 
                                              selected = 1)
                           )
                     
                 )
          )
        ),
        fluidRow(
          column(4,
                 box(title = "Normalisation Factors", width = 12, status = "info",
                    tableOutput("normFactorGenes")
                 )
                 ),
          column(4,
                 infoBoxOutput("Total_number_of_genes_HTSeq", width = 12)
          ),
          column(4,
                 infoBoxOutput("Number_of_genes_Filt_counts_Gene", width = 12)
          )
          
        ),
        br(),
        box(width="4", align="center", status="warning",
            actionButton("DEAnalysisOnlyGeneButt", "Differential Expression Analysis")
        ),
        br(),
        fluidRow(
              box(title = "Boxplot of raw counts", width = 12, status = "success",
                  highchartOutput("geneDensplot",height = "500px")
              ),
              box(title = "Boxplot of normalised counts", width = 12, status = "success",
                  highchartOutput("geneDensplotNormalised",height = "500px")
              ),
              box(title = "Differential expression analysis on genes", width = 12, status = "success",
                  DT::dataTableOutput("DEResultsTableOnlyOnGene")
              ),
              box(title = "Table for Up and Down regulated genes", width = 12, status = "success",
                  column(6,
                        highchartOutput("PlotUpDownForGene")
                  ),
                  column(6,
                        DT::dataTableOutput("TableUpDownForGene") 
                  )
              ),
              box(title = "Table for normalised counts", width = 12, status = "success",
                  DT::dataTableOutput("GeneNormalisedCounts")
              )
        )
        
)
