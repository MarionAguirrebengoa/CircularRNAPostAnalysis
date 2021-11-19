output$nbReplicats <- renderInfoBox({
  infoBox(
    "Replicats", length(unique(data.shiny$dataCircRNAPrediction$Replicat)), icon = icon("clone"),
    color = "blue"
  )
})

output$nbCdn <- renderInfoBox({
  infoBox(
    "Conditions", length(unique(data.shiny$dataCircRNAPrediction$Condition)), icon = icon("table"),
    color = "blue"
  )
})

output$nbSample <- renderInfoBox({
  infoBox(
    "Samples", length(unique(data.shiny$dataCircRNAPrediction$Sample)), icon = icon("male"),
    color = "blue"
  )
})
output$summary <- renderTable({
  withProgress(message = "Compute circular resume table", value = 0.1, {
    tryCatch({
    Sample <- unique(data.shiny$dataCircRNAPrediction$Sample)
    incProgress(0.3, detail = "Get Samples")
    # Car certaine jonction sont duppliqué de par leur annotation
    if (input$isStrandedData){
        DeleteDuplicateJunctions <- dcast(data.shiny$dataCircRNAPrediction, formula = Chromosome + Start +  End + Length + Strand + Sample + Condition + Replicat + NbSuppReads ~ .)
    }else{
        DeleteDuplicateJunctions <- dcast(data.shiny$dataCircRNAPrediction, formula = Chromosome + Start +  End + Length + Sample + Condition + Replicat + NbSuppReads ~ .)
    }
    dd = data.frame(Sample = Sample, Condition = sapply(X = Sample, function(x){return(unique(data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$Sample == x, ]$Condition))}),
                    Replicate = sapply(X = Sample, function(x){return(unique(data.shiny$dataCircRNAPrediction[data.shiny$dataCircRNAPrediction$Sample == x, ]$Replicat))}),
                    Observations = sapply(X = Sample, FUN = function(x){return(nrow(DeleteDuplicateJunctions[DeleteDuplicateJunctions$Sample == x,]))}),
                    NbTotalCircularRNAReads = sapply(X = Sample, FUN = function(x){return(sum(DeleteDuplicateJunctions[DeleteDuplicateJunctions$Sample == x,]$NbSuppReads))})
                    )
    incProgress(1, detail = "Build Table")
    return(dd)
    }, error = function(e){
        showNotification("Error in output$summary, can't access to data !", type = "error")
        return(NULL)
    })
   })
})

observeEvent(input$ExportToRDS, {
  saveRDS(data.shiny$dataCircRNAPrediction, "circRNAPrediction_save.rds")
})




output$summaryProcessing <- DT::renderDataTable(
    data.shiny$dataProcessingInfo[, c("NbReadsTotal", "NbReadsMappedBySTAR", "NbReadsMappedByBowtie2", "NbSplicedReadsSTAR", "NBSplicedReadsBowtie2", "NbChimericReads", "NbAnchors", "NbFilteredChimericReads", "NbUniqueChimericJunctions", "Sample", "Condition", "Replicat")],
    style="bootstrap",
    rownames = FALSE,
    escape = FALSE,
    extensions = 'Buttons',
    options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'))
)

output$MergingInformations <- DT::renderDataTable(
    data.shiny$dataProcessingInfo[, c("Sample", "Condition", "Replicat", "Mean", "Median", "Mode", "STDev", "PercentOfPairs")],
    style="bootstrap",
    rownames = FALSE,
    escape = FALSE,
    extensions = 'Buttons',
    options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'))
)

output$circRNATable = DT::renderDataTable(
  data.shiny$dataCircRNAPrediction,
  style="bootstrap",
  rownames = FALSE,
  escape = FALSE,
  extensions = 'Buttons',
  options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)

observeEvent(input$GenerateLinksGB, {
  withProgress(message = "Generate link to UCSC", value = 0.1, {
    generateLinkGB()
  incProgress(1, detail = "Done")
  })
})

generateLinkGB <- reactive({
  if(!is.null(data.shiny$dataCircRNAPrediction)){
    data.shiny$dataCircRNAPrediction$UCSC = paste0("<a href=\"http://genome.ucsc.edu/cgi-bin/hgTracks?org=", input$OrganismGB, "&db=", input$dataBaseGB, "&position=", data.shiny$dataCircRNAPrediction$Chromosome, ":", data.shiny$dataCircRNAPrediction$Start, "-", data.shiny$dataCircRNAPrediction$End, "\" target=\"_blank\">UCSC Genome Browser</a>")
  }
})


output$download_selected_rows <- downloadHandler(
  filename = "circRNA_predictions_selected_rows.csv",
  content = function(file){
    dataToDL = data.shiny$dataCircRNAPrediction[input$circRNATable_rows_selected, ]
    write.table(x = dataToDL, file, sep = ",", row.names = FALSE)
  }
)







