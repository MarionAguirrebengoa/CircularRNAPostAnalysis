output$Chromosomes <- renderUI({
  chromosomes <- as.vector(sort(unique(data.shiny$dataCircRNAPrediction$Chromosome)))
  selectInput(inputId = "ChromosomeList", label = "Chromosomes", choices = chromosomes, selectize = TRUE, multiple = TRUE)
})

filterDataCirc <- reactive({
    withProgress(message = "Compute circular resume table", value = 0.1, {
        incProgress(1/16, detail = "Chromosome filtering")
        if (!is.null(input$ChromosomeList)){
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$Chromosome %in% input$ChromosomeList, ]
        }
        incProgress(2/16, detail = "Replicate filtering")
        if (!is.null(input$ReplicateList)){
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$Replicat %in% input$ReplicateList, ]
        }
        incProgress(3/16, detail = "Condition filtering")
        if (!is.null(input$ConditionsList)){
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$Condition %in% input$ConditionsList, ]
        }
        incProgress(4/16, detail = "Sample filtering")
        if (!is.null(input$SamplesList)){
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$Sample %in% input$SamplesList, ]
        }
        incProgress(5/16, detail = "Length filtering")
        if (!is.na(input$minLength) && !is.na(input$maxLength)){
            data.shiny$dataCircRNAPrediction <- subset(data.shiny$dataCircRNAPrediction, data.shiny$dataCircRNAPrediction$Length >= input$minLength & data.shiny$dataCircRNAPrediction$Length <= input$maxLength)
        }
        incProgress(6/16, detail = "Signal filtering")
        if (!is.null(input$SignalList)){
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$Signal %in% input$SignalList, ]
        }
        incProgress(7/16, detail = "Supported reads filtering")
        if (!is.na(input$MinimalReadsSupported)){
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$NbSuppReads >= input$MinimalReadsSupported, ]                            
        }
        incProgress(8/16, detail = "Supported methods filtering")
        if (!is.na(input$MinimalSupportedMethods)){
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$supportedMethod >= input$MinimalSupportedMethods, ]
        }
        incProgress(9/16, detail = "Class filtering")
        if (!is.null(input$ClassesList)){
            print("Classes")
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$Class %in% input$ClassesList, ]
            }
        incProgress(10/16, detail = "Gene biotype filtering")
        if (!is.null(input$GeneBioTypeList)){
            print("GeneBioType")
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$GeneBioType %in% input$GeneBioTypeList, ]
        }
        incProgress(11/16, detail = "Gene symbol filtering")
        if (!is.null(input$GeneSymbolList)){
            print("GeneSymbol")
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$GeneSymbol %in% input$GeneSymbolList, ]
        }
        incProgress(12/16, detail = "Exact match filtering")
        if (input$ExactMatchList != "ANY"){
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$ExactMatch  == input$ExactMatchList, ]
        }
        incProgress(12/16, detail = "Chimeras filtering")
        if (input$isChimeraList != "ANY"){
            data.shiny$dataCircRNAPrediction <-  data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$isChimera  == input$isChimeraList, ]
        }
        incProgress(13/16, detail = "Gene ids filtering")
        if (!is.null(input$GeneIdList)){
            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$GeneId %in% input$GeneIdList, ]
        }
        incProgress(14/16, detail = "Use merged reads filtering")
        if (!is.null(input$useMerge) && input$useMerge != "ANY"){
            if(input$globalMerged){
                tableInter <- data.shiny$dataCircRNAPrediction %>% group_by(Chromosome, Start, End, Strand) %>% mutate(merged = any(merged))
                data.shiny$dataCircRNAPrediction <- as.data.frame(tableInter[tableInter$merged == input$useMerge, ])
            }else{
                data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$merged == input$useMerge, ]
            }
        }
        incProgress(15/16, detail = "Use unmerged reads filtering")
        if (!is.null(input$useUnmerge) && input$useUnmerge != "ANY"){
            if(input$globalUnmerged){
                tableInter <- data.shiny$dataCircRNAPrediction %>% group_by(Chromosome, Start, End, Strand) %>% mutate(unmerged = any(unmerged))
                data.shiny$dataCircRNAPrediction <- as.data.frame(tableInter[tableInter$unmerged == input$useUnmerge, ])
            }else{
                data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$unmerged == input$useUnmerge, ]
            }
        }
        incProgress(16/16, detail = "Done")
    })
})

observeEvent(input$filterButton, {
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    filterDataCirc()
  }
})

resetDataCirc <- reactive({
  data.shiny$dataCircRNAPrediction <- getDataCircFromInput()
})

observeEvent(input$ResetButton, {
  resetDataCirc()
})

output$Replicate <- renderUI({
  withProgress(message = "Getting replicates", value = 0.1, {
    Replicates <- as.vector(unique(data.shiny$dataCircRNAPrediction$Replicat))
    incProgress(1, detail = "Done")
  })
  selectInput(inputId = "ReplicateList", label = "Replicates", choices = Replicates, selectize = TRUE, multiple = TRUE)
})

output$Conditions <- renderUI({
  withProgress(message = "Getting conditions", value = 0.1, {
    Conditions <- as.vector(unique(data.shiny$dataCircRNAPrediction$Condition))
    incProgress(1, detail = "Done")
  })
  selectInput(inputId = "ConditionsList", label = "Conditions", choices = Conditions, selectize = TRUE, multiple = TRUE)

})

output$Samples <- renderUI({
  withProgress(message = "Getting samples", value = 0.1, {
    Samples <- as.vector(unique(data.shiny$dataCircRNAPrediction$Sample))
    incProgress(1, detail = "Done")
  })
  selectInput(inputId = "SamplesList", label = "Samples", choices = Samples, selectize = TRUE, multiple = TRUE)
})

output$LengthMin <- renderUI({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    withProgress(message = "Getting minimal length", value = 0.1, {
        minLength <- min(data.shiny$dataCircRNAPrediction$Length)
        incProgress(1, detail = "Done")
    })
    numericInput("minLength", label = paste("From, Current : ", minLength), value = NULL, min = minLength)
  }else{
    numericInput("minLength", label = "From", value = NULL)
  }
  })

output$LengthMax <- renderUI({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    withProgress(message = "Getting maximal length", value = 0.1, {
        maxLength <- max(data.shiny$dataCircRNAPrediction$Length)
        incProgress(1, detail = "Done")
    })
    numericInput("maxLength", label = paste("To, Current : ", maxLength), value = NULL, max = maxLength)
  }else{
    numericInput("maxLength", label = "To", value = NULL)
  }
})

output$Signal <- renderUI({
  withProgress(message = "Getting signals", value = 0.1, {
    Signals <- unique(data.shiny$dataCircRNAPrediction$Signal)
    incProgress(1, detail = "Done")
  })
    selectInput(inputId = "SignalList", label = "Signal", choices = Signals, selectize = TRUE, multiple = TRUE)
})

output$JunctionsReads <- renderUI({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
  withProgress(message = "Getting junction reads", value = 0.1, {
    minimalReads <- min(data.shiny$dataCircRNAPrediction$NbSuppReads[data.shiny$dataCircRNAPrediction$NbSuppReads != 0])
    incProgress(1, detail = "Done")
  })
    numericInput("MinimalReadsSupported", label = paste("Minimal Reads Supported, Current : ", minimalReads),  value = NULL, min = 0)
  }else{
    numericInput("MinimalReadsSupported", label = "Minimal Reads Supported (for at least one method)", value = NULL, min = 0)
  }
})

output$SupportedMethods <- renderUI({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    withProgress(message = "Getting supported methods", value = 0.1, {
        minimalMethods <- min(data.shiny$dataCircRNAPrediction$supportedMethod)
        incProgress(1, detail = "Done")
  })
    numericInput("MinimalSupportedMethods", label = paste("Minimal Supported Methods, Current : ", minimalMethods), value = NULL, min = 0, max = 3)
  }else{
    numericInput("MinimalSupportedMethods", label = "Minimal Supported Methods", value = NULL, min = 0, max = 3)
  }
})

output$Classes <- renderUI({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    withProgress(message = "Getting classes", value = 0.1, {
        ClassesList = unique(data.shiny$dataCircRNAPrediction$Class)
        incProgress(1, detail = "Done")
  })
    selectInput(inputId = "ClassesList", label = "Classes", choices = ClassesList, selectize = TRUE, multiple = TRUE)
  }else{
    selectInput(inputId = "ClassesList", label = "Classes", choices = NULL, selectize = TRUE, multiple = TRUE, width = '100%')
  }
})




output$GeneBiotypes <- renderUI({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    withProgress(message = "Getting gene biotypes", value = 0.1, {
        GeneBioTypeSet <- unique(data.shiny$dataCircRNAPrediction$GeneBioType) %>% as.character()
        GeneBioTypeSet <- GeneBioTypeSet[!is.na(GeneBioTypeSet)]
        incProgress(1, detail = "Done")
  })
    selectInput(inputId = "GeneBioTypeList", label = "GeneBiotypes", choices = GeneBioTypeSet, selectize = TRUE, multiple = TRUE, width = '100%')
  }else{
    selectInput(inputId = "GeneBioTypeList", label = "GeneBiotypes", choices = NULL, selectize = TRUE, multiple = TRUE, width = '100%')
  }
})



 output$Genes <- renderUI({
   if (!is.null(data.shiny$dataCircRNAPrediction)){
    withProgress(message = "Getting genes symbol", value = 0.1, {
        GeneSymbolList = unique(data.shiny$dataCircRNAPrediction$GeneSymbol) %>% as.character()
        GeneSymbolList <- GeneSymbolList[!is.na(GeneSymbolList)]
        incProgress(1, detail = "Done")
  })
    selectInput(inputId = "GeneSymbolList", label = "GeneSymbol", choices = GeneSymbolList, selectize = TRUE, multiple = TRUE)
   }else{
     selectInput(inputId = "GeneSymbolList", label = "GeneSymbol", choices = NULL, selectize = TRUE, multiple = TRUE)
   }

 })

 output$ExactMatch <- renderUI({
   if (!is.null(data.shiny$dataCircRNAPrediction)){
        selectInput(inputId = "ExactMatchList", label = "ExactMatch", choices = c(TRUE, FALSE, "ANY"), selectize = TRUE, multiple = FALSE, selected = "ANY")
   }
 })

 output$isChimera <- renderUI({
   if (!is.null(data.shiny$dataCircRNAPrediction)){
        selectInput(inputId = "isChimeraList", label = "isChimera", choices = c(TRUE, FALSE, "ANY"), selectize = TRUE, multiple = FALSE, selected = "ANY")
   }
 })

 output$useMerge <- renderUI({
     if (!is.null(data.shiny$dataCircRNAPrediction) && ("merged" %in% colnames(data.shiny$dataCircRNAPrediction))){
        selectInput(inputId = "useMerge", label = "use merged reads", choices = c(TRUE, FALSE, "ANY"), selectize = TRUE, multiple = FALSE, selected = "ANY")
     }else{
        selectInput(inputId = "useMerge", label = "use merged reads", choices = NULL, selectize = TRUE, multiple = FALSE)
   }
 })

 output$useUnmerge <- renderUI({
     if (!is.null(data.shiny$dataCircRNAPrediction) && ("unmerged" %in% colnames(data.shiny$dataCircRNAPrediction))){
        selectInput(inputId = "useUnmerge", label = "use unmerged reads", choices = c(TRUE, FALSE, "ANY"), selectize = TRUE, multiple = FALSE, selected = "ANY")
     }else{
         selectInput(inputId = "useUnmerge", label = "use unmerged reads", choices = NULL, selectize = TRUE, multiple = FALSE)
     }
 })


 output$GeneId <- renderUI({
   if (!is.null(data.shiny$dataCircRNAPrediction)){
     withProgress(message = "Getting genes ids", value = 0.1, {
        GeneIdList <- unique(data.shiny$dataCircRNAPrediction$GeneId)
        GeneIdList <- GeneIdList[!is.na(GeneIdList)]
        incProgress(1, detail = "Done")
  })
    selectInput(inputId = "GeneIdList", label = "GeneId", choices = GeneIdList, selectize = TRUE, multiple = TRUE)
   }else{
     selectInput(inputId = "GeneIdList", label = "GeneId", choices = NULL, selectize = TRUE, multiple = TRUE)
 }
 })

# For repeat regions


  observeEvent(input$filechoose, {
    data.shiny$pathToRepeatRegionsFile <- file.choose()
  })

   output$filechosen <- renderText({
      if (is.null(data.shiny$pathToRepeatRegionsFile)){
        "Nothing selected"
      }else{
        data.shiny$pathToRepeatRegionsFile
      }
   })

   observeEvent(input$removeRepeatRegions, {
        if (!is.null(data.shiny$pathToRepeatRegionsFile)){
            if (file.exists(data.shiny$pathToRepeatRegionsFile)){
                withProgress(message = "Import File", value = 0.1, {
                    tryCatch({
                        incProgress(0.2, detail = "Converting to GRanges")
                        repeatRegion <- import.bed(data.shiny$pathToRepeatRegionsFile)
                        incProgress(0.4, detail = "Converting circular coordinates to GRanges")
                        startCoordinates <- GRanges(seqnames = data.shiny$dataCircRNAPrediction$Chromosome,
                                                    ranges = IRanges(data.shiny$dataCircRNAPrediction$Start + 1, data.shiny$dataCircRNAPrediction$Start + 2),
                                                    strand = data.shiny$dataCircRNAPrediction$Strand)
                        endCoordinates <- GRanges(seqnames = data.shiny$dataCircRNAPrediction$Chromosome,
                                                    ranges = IRanges(data.shiny$dataCircRNAPrediction$End-1, data.shiny$dataCircRNAPrediction$End),
                                                    strand = data.shiny$dataCircRNAPrediction$Strand)

                        startOverlapRepeatRegions <- queryHits(findOverlaps(query = startCoordinates, subject = repeatRegion, select = "all"))
                        endOverlapRepeatRegions <- queryHits(findOverlaps(query = endCoordinates, subject = repeatRegion, select = "all"))
                        incProgress(0.8, detail = "Merging")
                        allOverlapRepeatRegions <- union(startOverlapRepeatRegions, endOverlapRepeatRegions)
                        incProgress(0.9, detail = "Remove elements")
                        if (length(allOverlapRepeatRegions) > 0){
                            data.shiny$dataCircRNAPrediction <- data.shiny$dataCircRNAPrediction[-allOverlapRepeatRegions, ]
                        }
                        incProgress(1, detail = "Done")
                    }, error = function(e){
                        showNotification("Error removing circular matching repeat regions !", type = "error")
                    })
                })
            }else{
                showNotification("Repeat file does not exist !", type = "error")
            }
        }
   })