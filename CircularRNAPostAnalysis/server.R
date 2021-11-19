# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#    http://shiny.rstudio.com/
#
####Library########
library(shiny)
library(dplyr)
library(DT)
library(ggplot2)
library(highcharter)
library(DESeq2)
library(reshape2)
library(biomaRt)
library(rtracklayer)
library(GenomicRanges)
###################

options(shiny.maxRequestSize=50*1024^2)
options(shiny.sanitize.errors = FALSE)

ensembl_mart <- tryCatch({
  useMart("ensembl")
}, error=function(err){
  return(NULL)
})

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  data.shiny <- reactiveValues(dataCircRNAPrediction = NULL,
                               dataProcessingInfo = NULL,
                               dataGeneCount = NULL,
                               countMatrix = NULL,
                               coldata = NULL,
                               TotalNumberOfGenesBeforeFiltering = NULL,
                               TotalNumberOfGenesFromHTSeqCounts = NULL,
                               circRNANormalisedCounts = data.frame(),
                               geneNormalisedCounts = data.frame(),
                               circRNAFinalCounts = NULL,
                               DEResults = NULL,
                               
                               dataForKSTest = NULL,
                               KSTestResults = NULL,
                               circRNAGloablNormalisedCounts = NULL,
                               DEResultsOnlyOnCirc = NULL,
                               DEResultsOnlyOnGene = NULL,
                               FilteredJunctionsFromCirc = data.frame(),
                               ensembl = ensembl_mart,
                               UpDownCirc = data.frame(),
                               UpDownGene = data.frame(),
                               UpDownJunction = data.frame(),
                               
                               dataForClust = NULL,
                               CAHResult = NULL,
                               ResultatsFromCAH_points = NULL,
                               ResultatsFromCAH_Barycentres = NULL,
                               
                               DEResultsOnlyCircularJunctions = NULL,
                               dataUpAndDownForCirc = NULL,
                               dataUpAndDownForJunctions = NULL,
                               FilteredJunctionsFromCircDescriptions = NULL,
                               JunctionNormalisedCounts = data.frame(),
                               
                               normFactorForGenes = NULL,

                               pathToRepeatRegionsFile = NULL
                               
  )

  source(file.path("server", "load_dataServer.R"), local = TRUE)$value
  source(file.path("server", "SummaryServer.R"), local = TRUE)$value
  source(file.path("server", "filterServer.R"), local = TRUE)$value
  source(file.path("server", "PrimaryStatServer.R"), local = TRUE)$value
  source(file.path("server", "DEAnalysisServerGene.R"), local = TRUE)$value
  source(file.path("server", "DEAnalysisServerCirc.R"), local = TRUE)$value
  source(file.path("server", "DEAnalysisServerCircGene.R"), local = TRUE)$value
  
  
  output$TestIfDeIsDoneOnGeneANDCircular <- reactive({

    return((!is.null(data.shiny$DEResultsOnlyOnGene) && !is.null(data.shiny$DEResultsOnlyOnCirc)))
  })

  outputOptions(output, "TestIfDeIsDoneOnGeneANDCircular", suspendWhenHidden = FALSE)
  })

