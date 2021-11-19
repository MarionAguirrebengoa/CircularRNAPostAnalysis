
computeCAH <- function(dataOnCirculars, dataOnGenes){
  if (!is.null(dataOnCirculars) && !is.null(dataOnGenes)){
    withProgress(message = "Compute CAH", value = 0.1, {
        tryCatch({
                incProgress(0.2, detail = "Get circular and genes expression")
                # On ne selctionne que ceux qui sont vue partouts
                dataGene <- dataOnGenes[rownames(dataOnGenes) %in% rownames(dataOnCirculars), ]
                dataCirc <- dataOnCirculars[rownames(dataOnCirculars) %in% rownames(dataGene), ]
                # On les mets dans le même ordre
                dataGene <- dataGene[rownames(dataCirc), ]
                dataForClust <- data.frame(GeneId = rownames(dataCirc), Log2FCCirc = dataCirc$log2FoldChange, Log2FCGene = dataGene$log2FoldChange)

                incProgress(0.5, detail = "Compute clusters")
                distances <- dist(dataForClust[ ,c("Log2FCCirc", "Log2FCGene")], method = "euclidean")
                cah.test <- hclust(distances, method = "ward.D2")
                data.shiny$dataForClust <- dataForClust

                incProgress(1, detail = "Done")
                return(cah.test)
        }, error = function(e){
            showNotification("Error in getTableOfCommon!", type = "error")
            data.shiny$dataForClust <- NULL
            return(NULL)
        })
    })
  }
}


launchCAH <- observeEvent(input$CAH_Button, {
  data.shiny$CAHResult = computeCAH(data.shiny$DEResultsOnlyOnCirc, data.shiny$DEResultsOnlyOnGene);
})


output$Dendrogramme_CAH <- renderPlot({
  if (!is.null(data.shiny$CAHResult)){
    plot(data.shiny$CAHResult, labels = FALSE, main = "Dendrogramme CAH", xlab = "", ylab = "", sub = "", axes = FALSE, hang = -1)
  }
})

output$Intertie <- renderPlot({
  if (!is.null(data.shiny$CAHResult)){
    inertie <- sort(data.shiny$CAHResult$height, decreasing = TRUE)
    plot(inertie[1:20], type = "s", xlab = "Nombre de classes", ylab = "Inertie")
  }
})


ReloadPlotOfCAH <- observeEvent(input$CAH_ScatterPlot_Button, {
    if (!is.null(data.shiny$CAHResult) && !is.null(data.shiny$dataForClust)){
        withProgress(message = "Display Scatter plot", value = 0.1, {
            incProgress(0.2, detail = "Get classes")
            groups <- cutree(data.shiny$CAHResult, k = input$NbClasses)
            dataForPlot <- data.shiny$dataForClust
            dataForPlot$Class <- groups
            data.shiny$ResultatsFromCAH_points <- dataForPlot

            incProgress(0.6, detail = "Get barycentres")
            Barycentres <- aggregate(dataForPlot[,c("Log2FCCirc", "Log2FCGene")], by=list(dataForPlot$Class), FUN = mean)
            colnames(Barycentres) <- c("Class", "Log2FCCirc", "Log2FCGene")
            data.shiny$ResultatsFromCAH_Barycentres <- Barycentres
            incProgress(1, detail = "Done")
        })
    }
})
 
 output$CAH_ScatterPlot <- renderHighchart({
   if (!is.null(data.shiny$ResultatsFromCAH_points) && !is.null(data.shiny$ResultatsFromCAH_Barycentres)){
        dataForClust <- data.shiny$ResultatsFromCAH_points %>% as.data.frame()
        Barycentres <- data.shiny$ResultatsFromCAH_Barycentres %>% as.data.frame()

        selectedGenes <- data.shiny$ResultatsFromCAH_points[data.shiny$ResultatsFromCAH_points$GeneId %in% input$SelectGeneId, ]

        x <- c("GeneId : ", "Log2FCCirc : ", "Log2FCGene : ")
        y <- c(sprintf("{point.%s}", c("GeneId")), sprintf("{point.%s:.2f}", c("Log2FCCirc", "Log2FCGene")))
        tltip <- tooltip_table(x, y)

        hc <- hchart(dataForClust, "scatter", hcaes(x = Log2FCCirc, y = Log2FCGene, group = Class)) %>%
        hc_tooltip(useHTML = TRUE, headerFormat = "", pointFormat = tltip) %>%
        hc_add_series(Barycentres, name="Barycentres", hcaes(x = Log2FCCirc, y = Log2FCGene), type = "scatter", marker = list(symbol = "cross", enabled = TRUE)) %>%
        hc_exporting(enabled = TRUE, filename = "CAH_ScatterPlot") %>%
        hc_chart(zoomType = "xy")
        if(length(input$SelectGeneId) > 0){
            hc <- hc %>% hc_add_series(selectedGenes, names="Selected_Genes", hcaes(x = Log2FCCirc, y = Log2FCGene), type = "scatter", marker = list(symbol = "triangle", enabled = TRUE), radius = 9)
        }
        hc
   }
 })

 output$download_scatter_CAH_points <- downloadHandler(
   filename = "CAH_results_by_class.csv",
   content = function(file){
     write.table(x = data.shiny$ResultatsFromCAH_points, file, sep = ",")
   }
 )

 output$download_scatter_CAH_barycentres <- downloadHandler(
   filename = "CAH_results_Barycentres_by_class.csv",
   content = function(file){
     write.table(x = data.shiny$ResultatsFromCAH_Barycentres, file, sep = ",")
   }
 )

output$SelectGeneId <- renderUI({
  GeneIdList = data.shiny$ResultatsFromCAH_points$GeneId
  selectInput(inputId = "SelectGeneId", label = "Select genes", choices = GeneIdList, selectize = TRUE, multiple = TRUE)
})
