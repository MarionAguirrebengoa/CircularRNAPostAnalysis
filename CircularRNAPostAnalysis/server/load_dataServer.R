

# Permet à partir de du fichier RDS, d'importer les données : On met les Na à 0, on construit la colonne NbSupported, et on met les mauvaus exacts Match (issu des non-annotated car pas dans le GTF: chrM rDNA) à FALSE
getDataCircFromInput <- reactive({
  if (!is.null(input$circRNAPred$datapath)){
    tryCatch({
        withProgress(message = "Import circular RNA predictions", value = 0.1, {
            data.shiny$dataCircRNAPrediction <- readRDS(input$circRNAPred$datapath)
            incProgress(0.4, detail = "Build Table")
            data.shiny$dataCircRNAPrediction$Junctions_reads_CircExplorer2[is.na(data.shiny$dataCircRNAPrediction$Junctions_reads_CircExplorer2)] <-0
            data.shiny$dataCircRNAPrediction$Junctions_reads_FindCirc[is.na(data.shiny$dataCircRNAPrediction$Junctions_reads_FindCirc)] <-0
            data.shiny$dataCircRNAPrediction$Junctions_reads_circRNA_finder[is.na(data.shiny$dataCircRNAPrediction$Junctions_reads_circRNA_finder)] <-0
            data.shiny$dataCircRNAPrediction <-  data.shiny$dataCircRNAPrediction %>% rowwise() %>% mutate(NbSuppReads = max(Junctions_reads_FindCirc, Junctions_reads_circRNA_finder, Junctions_reads_CircExplorer2)) %>% as.data.frame()
            # On s'assure qu'il s'agisse bien de character pour les geneIds
            data.shiny$dataCircRNAPrediction$GeneId <- as.character(data.shiny$dataCircRNAPrediction$GeneId)
            # On crée directement coldata
            incProgress(0.8, detail = "Creating coldata")
            data.shiny$coldata <- subset.data.frame(data.shiny$dataCircRNAPrediction, select = c("Sample", "Condition", "Replicat")) %>% distinct(Sample, Condition, Replicat) %>% as.data.frame()
            if (!is.null(data.shiny$dataGeneCount)){
                data.shiny$dataGeneCount <- data.shiny$dataGeneCount[, c(data.shiny$coldata$Sample) ]
            }
            incProgress(1, detail = "Done")
        })
    }, error = function(e){
      showNotification("Error in getDataCircFromInput, can't load data !", type = "error")
      data.shiny$dataCircRNAPrediction <- NULL
      data.shiny$coldata <- NULL
    })
    return(data.shiny$dataCircRNAPrediction)
  }else{
    showNotification("you don't load circular RNA predictions data !", type = "message")
  }
})

# Permet de récupérer les process
getDataProcessFromInput <- reactive({
  if (!is.null(input$ProcessInfo$datapath)){
    tryCatch({
      withProgress(message = "Import processing information", value = 0.1, {
        data.shiny$dataProcessingInfo <- readRDS(input$ProcessInfo$datapath)
        incProgress(1, detail = "Done")
      })
      }, error = function(e){
      showNotification("Error in getDataProcessFromInput, can't load data !", type = "error")
        data.shiny$dataProcessingInfo <- NULL
    })
    return(data.shiny$dataProcessingInfo)
  }else{
    showNotification("you don't load process informations data !", type = "message")
  }
})

# Pour récupérer les Reads conts à partir du fichier renvoyé du pipeline issu de HtSeq-count
getGeneReadCountFromInput <- reactive({
  if (!is.null(input$GeneReadCount$datapath)){
    tryCatch({
        withProgress(message = "Import gene counts information", value = 0.1, {
            dataGene <- readRDS(input$GeneReadCount$datapath)
            # On transforme en matrice avec les GeneIds en rownames
            incProgress(0.6, detail = "Format data")
            data.shiny$dataGeneCount <- data.frame(dataGene[, -1], row.names = as.character(dataGene[, 1]))
            if (!is.null(data.shiny$coldata)){
                data.shiny$dataGeneCount <- data.shiny$dataGeneCount[, c(data.shiny$coldata$Sample) ]
            }
            incProgress(1, detail = "Get total number of genes")
            data.shiny$TotalNumberOfGenesFromHTSeqCounts <- nrow(data.shiny$dataGeneCount)
        })
    }, error = function(e){
      showNotification("Error in getGeneReadCountFromInput, can't access to data !", type = "error")
      data.shiny$dataGeneCount <- NULL
    })
    return(data.shiny$dataGeneCount)
  }else{
    showNotification("you don't load process gene counts data !", type = "message")
  }
})

getColdataFromInput <- reactive({
        tryCatch({
              withProgress(message = "Import coldata information from phenotype", value = 0.1, {
                data.shiny$coldata <- read.table(file = input$PhenotypeInfo$datapath, header = TRUE, sep = ",")
                incProgress(1, detail = "Done")
              })
        }, error = function(e){
              showNotification("Error in getColdataFromInput, can't access to data !", type = "error")
              data.shiny$coldata <- NULL
        })
  
})


observeEvent(input$PhenotypeInfo$datapath, {
  getColdataFromInput()
})

observeEvent(input$GeneReadCount$datapath, {
  getGeneReadCountFromInput()
})

observeEvent(input$ProcessInfo$datapath, {
  getDataProcessFromInput()
})

observeEvent(input$circRNAPred$datapath, {
  getDataCircFromInput()
})
