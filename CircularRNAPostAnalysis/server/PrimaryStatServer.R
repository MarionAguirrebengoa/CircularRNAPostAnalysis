
createtabForPlot <- eventReactive(input$PrintPlot, {
    withProgress(message = "Get classes", value = 0.1, {
        dataForClassPlot <- data.shiny$dataCircRNAPrediction %>% group_by(Replicat,Condition,Class) %>% summarise(n = n()) %>% mutate(freq = (n/sum(n))*100) 
    incProgress(1, detail = "Done")
    })
    return(dataForClassPlot)
})




output$ClassPlot <- renderHighchart({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    myClassPlotData = createtabForPlot() %>% mutate(Series = paste(Condition, "-", Replicat))
    x <- c("Replicat : ", "Condition : ", "Frequency : ", "Raw value : ")
    y <- c( sprintf("{point.%s}", c("Replicat", "Condition")), sprintf("{point.%s:.2f}", c("freq", "n")))
    tltip <- tooltip_table(x, y)
    hchart(myClassPlotData, "column", hcaes(x = Class, y = freq, group = Series)) %>% hc_tooltip(useHTML = TRUE, headerFormat = "", pointFormat = tltip) %>%  
    hc_title(text = "Classes Distribution") %>% 
    hc_yAxis(title = list(text = "Frequency"), labels = list(format = "{value}%")) %>%
    hc_add_theme(hc_theme_smpl()) %>%
    hc_exporting(enabled = TRUE, filename = "Class_Distribution_Plot", formAttributes = list(target = "_blank"))  %>%
    hc_chart(zoomType = "xy")
  }
})



output$ClassTable <- DT::renderDataTable(
  createtabForPlot(),
  style="bootstrap",
  rownames = FALSE,
  escape = FALSE,
  extensions = 'Buttons',
  options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)


getLengthTable <- function(data){ 
    startVector = c(0,100,250,500,1000,5000,10000,50000, 100000)
    endVector = c(100, 250,500,1000,5000,10000,50000,100000,200000)
    seqVec <- Vectorize(seq.default, vectorize.args = c("from", "to"))
    Cath = seqVec(from = startVector, to = endVector-1, by = 1)
    table = lapply(Cath, FUN = function(x){return(nrow(data[data$Length %in% x, ]))})
    superiorToMax = nrow(data[data$Length >= tail(endVector, n=1), ])
    res = data.frame(Intervals = c(paste0("[", startVector,"-", endVector, "]"), paste0("[", tail(endVector, n=1),"-[")), Values = c(as.vector(unlist(table)), superiorToMax))
    res$Pourcent = prop.table(res$Values)*100
    res$Replicat = rep(data[1,]$Replicat, length(startVector) + 1) 
    res$Condition = rep(data[1,]$Condition, length(startVector) + 1)
  return(res)
}


createtabForPlotLength <- eventReactive(input$PrintPlotLength, {
  withProgress(message = "Get lengths", value = 0.1, {
    DeleteDuplicateJunctions <- dcast(data.shiny$dataCircRNAPrediction, formula = Chromosome + Start +  End + Length + Strand + Sample + Condition + Replicat ~ .)
    firstStep <- split(DeleteDuplicateJunctions, list(DeleteDuplicateJunctions$Replicat, DeleteDuplicateJunctions$Condition))
    incProgress(0.5, detail = "Preparing data")
    dataForLengthPlot <- lapply(firstStep, FUN = function(x){ return(getLengthTable(x))}) %>% bind_rows()
    incProgress(1, detail = "Done")
})
  return(dataForLengthPlot)
  
})



output$LengthPlot <- renderHighchart({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    myLengthPlotData <- createtabForPlotLength() %>% mutate(Series = paste(Condition, "-", Replicat))
    x <- c("Replicat : ", "Condition : ", "Percentage : ", "Raw value : ")
    y <- c( sprintf("{point.%s}", c("Replicat", "Condition")), sprintf("{point.%s:.2f}", c("Pourcent", "Values")))
    tltip <- tooltip_table(x, y)
    hchart(myLengthPlotData, "column", hcaes(x = Intervals, y = Pourcent, group = Series)) %>% hc_tooltip(useHTML = TRUE, headerFormat = "", pointFormat = tltip) %>%  
        hc_title(text = "Length Distribution") %>% 
        hc_yAxis(title = list(text = "Poucentage"), labels = list(format = "{value}%")) %>%
        hc_add_theme(hc_theme_smpl()) %>%
        hc_exporting(enabled = TRUE, filename = "Length_Distribution_Plot", formAttributes = list(target = "_blank")) %>%
        hc_chart(zoomType = "xy")
 }
})


output$LengthTable = DT::renderDataTable(
  createtabForPlotLength(),
  style="bootstrap",
  rownames = FALSE,
  escape = FALSE,
  extensions = 'Buttons',
  options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)



createTabForSignal <- eventReactive(input$PrintPlotSignal, {
  withProgress(message = "Get signals", value = 0.1, {
    DeleteDuplicateJunctions <- dcast(data.shiny$dataCircRNAPrediction, formula = Chromosome + Start +  End + Length + Strand + Sample + Condition + Replicat + Signal ~ .)
    incProgress(0.5, detail = "...")
    dataForSignal <- DeleteDuplicateJunctions %>% group_by(Replicat, Condition, Signal) %>% summarise(n=n()) %>% mutate(freq = (n/sum(n))*100)
    incProgress(1, detail = "Done")
  })
  return(dataForSignal[order(dataForSignal$freq, decreasing = TRUE), ])
})


output$SignalPlot <- renderHighchart({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    mySignalPlotData <- createTabForSignal() %>% mutate(Series = paste(Condition, "-", Replicat))
    x <- c("Replicat : ", "Condition : ", "Frequency : ", "Raw value : ")
    y <- c( sprintf("{point.%s}", c("Replicat", "Condition")), sprintf("{point.%s:.2f}", c("freq", "n")))
    tltip <- tooltip_table(x, y)
    hchart(mySignalPlotData, "column", hcaes(x = Signal, y = freq, group = Series)) %>% hc_tooltip(useHTML = TRUE, headerFormat = "", pointFormat = tltip) %>%  
        hc_title(text = "Signal Distribution") %>% 
        hc_yAxis(title = list(text = "Frequency"), labels = list(format = "{value}%")) %>%
        hc_add_theme(hc_theme_smpl()) %>%
        hc_exporting(enabled = TRUE, filename = "Signal_Distribution_Plot", formAttributes = list(target = "_blank")) %>%
        hc_chart(zoomType = "xy")
  }
})

output$SignalTable <- DT::renderDataTable(
  createTabForSignal(),
  style="bootstrap",
  rownames = FALSE,
  escape = FALSE,
  extensions = 'Buttons',
  options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)


createTableForSuppReads <- eventReactive(input$PrintPlotSuppReads, {
  withProgress(message = "Get supported reads", value = 0.1, {
        DeleteDuplicateJunctions = dcast(data.shiny$dataCircRNAPrediction, formula = Chromosome + Start +  End + Length + Strand + Sample + Condition + Replicat + NbSuppReads ~ .)
        dataForSuppReads = DeleteDuplicateJunctions %>% group_by(Replicat, Condition, NbSuppReads) %>% summarise(n = n()) %>% mutate(freq = (n/sum(n))*100)
        incProgress(1, detail = "Done")
  })
  return(dataForSuppReads[order(dataForSuppReads$freq, decreasing = TRUE), ])
})

output$SuppReadsPlot <- renderHighchart({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    mySuppReadsPlotData = createTableForSuppReads() %>% mutate(Series = paste(Condition, "-", Replicat))
    x <- c("Replicat : ", "Condition : ", "Frequency : ", "Raw value : ")
    y <- c(sprintf("{point.%s}", c("Replicat", "Condition")), sprintf("{point.%s:.2f}", c("freq", "n")))
    tltip <- tooltip_table(x, y)
    hchart(mySuppReadsPlotData, "column", hcaes(x = NbSuppReads, y = freq, group = Series)) %>% hc_tooltip(useHTML = TRUE, headerFormat = "", pointFormat = tltip) %>%  
        hc_title(text = "Supported Reads Distribution") %>% 
        hc_yAxis(title = list(text = "Frequency"), labels = list(format = "{value}%")) %>%
        hc_add_theme(hc_theme_smpl()) %>%
        hc_exporting(enabled = TRUE, filename = "Supported_Reads_Distribution_Plot", formAttributes = list(target = "_blank")) %>%
        hc_chart(zoomType = "xy")
  }
})


output$SuppReadsTable <- DT::renderDataTable(
  createTableForSuppReads(),
  style="bootstrap",
  rownames = FALSE,
  escape = FALSE,
  extensions = 'Buttons',
  options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)



createTabForSuppMeth <- eventReactive(input$PrintPlotSuppMeth, {
  withProgress(message = "Compute supported methods", value = 0.1, {
    DeleteDuplicateJunctions = dcast(data.shiny$dataCircRNAPrediction, formula = Chromosome + Start +  End + Length + Strand + Sample + Condition + Replicat + supportedMethod ~ .)
    dataForSuppMeth = DeleteDuplicateJunctions %>% group_by(Replicat, Condition, supportedMethod) %>% summarise(n = n()) %>% mutate(freq = (n/sum(n))*100)
    incProgress(1, detail = "Done")
 })
  return(dataForSuppMeth)
})


output$SuppMethPlot <- renderHighchart({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    mySuppMethPlotData <- createTabForSuppMeth() %>% mutate(Series = paste(Condition, "-", Replicat))
    x <- c("Replicat : ", "Condition : ", "Frequency : ", "Raw value : ")
    y <- c( sprintf("{point.%s}", c("Replicat", "Condition")), sprintf("{point.%s:.2f}", c("freq", "n")))
    tltip <- tooltip_table(x, y)
    hchart(mySuppMethPlotData, "column", hcaes(x = supportedMethod, y = freq, group = Series)) %>% hc_tooltip(useHTML = TRUE, headerFormat = "", pointFormat = tltip) %>%  
        hc_title(text = "Supported Methods Distribution") %>% 
        hc_yAxis(title = list(text = "Frequency"), labels = list(format = "{value}%")) %>%
        hc_add_theme(hc_theme_smpl()) %>%
        hc_exporting(enabled = TRUE, filename = "Supported_Meth_Distribution_Plot", formAttributes = list(target = "_blank")) %>%
        hc_chart(zoomType = "xy")
  }
})

output$SuppMethTable <- DT::renderDataTable(
  createTabForSuppMeth(),
  style="bootstrap",
  rownames = FALSE,
  escape = FALSE,
  extensions = 'Buttons',
  options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)


createTableForToolsStatistics <- eventReactive(input$PrintPlotToolsStatistics, {
  withProgress(message = "get tools statistics", value = 0.1, {
    DeleteDuplicateJunctions <- dcast(data.shiny$dataCircRNAPrediction, formula = Chromosome + Start +  End + Length + Strand + Sample + Junctions_reads_FindCirc + Junctions_reads_circRNA_finder + Junctions_reads_CircExplorer2 + Condition + Replicat ~ .)
    DataForToolsStats <- DeleteDuplicateJunctions %>% group_by(Replicat, Condition) %>% summarise(Find_circ = sum(Junctions_reads_FindCirc), CircRNA_finder = sum(Junctions_reads_circRNA_finder), CircExplorer2 = sum(Junctions_reads_CircExplorer2))
    incProgress(1, detail = "Done")
    })
    return(DataForToolsStats)
})


output$ToolsStatisticsPlot <- renderHighchart({
  if (!is.null(data.shiny$dataCircRNAPrediction)){
    myDataForToolsStatsPlot <- createTableForToolsStatistics() %>% mutate(Series = paste(Condition, "-", Replicat)) %>% melt(id.vars = c("Series"), variable.name = "Tool", measure.vars = c("Find_circ", "CircRNA_finder", "CircExplorer2"), value.name = "Reads")
    hchart(myDataForToolsStatsPlot, "column", hcaes(x = Tool, y = Reads, group = Series)) %>%  
        hc_title(text = "Tools detection statistics") %>% 
        hc_yAxis(title = list(text = "Number of detected reads")) %>%
        hc_add_theme(hc_theme_smpl()) %>%
        hc_exporting(enabled = TRUE, filename = "Tools_Statistics_Plot", formAttributes = list(target = "_blank")) %>%
        hc_chart(zoomType = "xy")
  }
})






