
output$mainEffectCirc <- renderUI({
  selectInput(inputId = "mainEffectCirc", label = "main effect to test", choices = colnames(data.shiny$coldata), selectize = TRUE, multiple = FALSE)
})

output$batchEffectCirc <- renderUI({
  selectInput(inputId = "batchEffectCirc", label = "batch effect", choices = colnames(data.shiny$coldata), selectize = TRUE, multiple = TRUE)
})

output$Condition_1_circ <- renderUI({
  if(!is.null(input$mainEffectCirc)){
    selectInput(inputId = "Condition_1_circ", label = "Condition 1", choices = data.shiny$coldata[[input$mainEffectCirc]], selectize = TRUE, multiple = FALSE)
  }
})

output$Condition_2_circ <- renderUI({
  if(!is.null(input$mainEffectCirc)){
    selectInput(inputId = "Condition_2_circ", label = "Condition 2", choices = data.shiny$coldata[[input$mainEffectCirc]], selectize = TRUE, multiple = FALSE)
  }
})

NormalisationCirc <- function(CountMatrix, coldata){
      dds <- NULL
      withProgress(message = "Normalise counts", value = 0.1, {
        tryCatch({
                incProgress(0.2, detail = "Get design to prepare object")
                Design <- paste("~", paste(paste(input$mainEffectCirc, collapse = " + "), paste(input$batchEffectCirc, collapse = " + "), sep = ifelse(!is.null(input$batchEffectCirc), yes = " + ", no = "")))
                # factor from genes
                if (input$NormMethod_Circ == 1){
                        incProgress(0.4, detail = "From genes")
                        # On prépare l'objet DeSeq
                        dds = DESeqDataSetFromMatrix(countData = CountMatrix, colData = coldata, design = formula(Design))
                        # On fait la normalisation avec le facteur utilisé pour les gènes
                        sizeFactors(dds) = data.shiny$normFactorForGenes
                }
                # CPM normalisation
                else if(input$NormMethod_Circ == 2){
                    incProgress(0.4, detail = "CPM")
                    NbMappedReads = colSums(CountMatrix)
                    CPMFactor = NbMappedReads/input$normFactCPM_Circ
                    # On prépare l'objet DeSeq
                    dds = DESeqDataSetFromMatrix(countData = CountMatrix, colData = coldata, design = formula(Design))
                    # On met le sizeFactor en CPM
                    sizeFactors(dds) = CPMFactor
                }else if(input$NormMethod_Circ == 3){
                    # On fait un objet classique
                    incProgress(0.4, detail = "Deseq")
                    dds = DESeqDataSetFromMatrix(countData = CountMatrix, colData = coldata, design = formula(Design))
                    dds <- estimateSizeFactors(dds)
                }else{
                    incProgress(0.4, detail = "RPM")
                    NbMappedReads = data.shiny$dataProcessingInfo[match(coldata$Sample, data.shiny$dataProcessingInfo$Sample),]$NbReadsMappedBySTAR
                    CPMFactor = NbMappedReads/input$normFactCPM_Circ
                    dds = DESeqDataSetFromMatrix(countData = CountMatrix, colData = coldata, design = formula(Design))
                    sizeFactors(dds) = CPMFactor
                }
            incProgress(1, detail = "Normalisation done")
        }, error = function(e){
        showNotification("Error in NormalisationCirc, can't normalise data !", type = "error")
        showNotification(e)
        dds <- NULL
    })
  })
    return(dds)
}

observeEvent(input$NormMethod_Circ, {
  if (input$NormMethod_Circ == 1 && is.null(data.shiny$normFactorForGenes)){
    shinyjs::disable("DEAnalysisOnlyCircButt")
    shinyjs::disable("DEAnalysisOnlyJunctionsButt")
  }else{
    shinyjs::enable("DEAnalysisOnlyCircButt")
    shinyjs::enable("DEAnalysisOnlyJunctionsButt")
  }
  
})



# Permet de créer la matrice de comptage
# On récupère les colonnes Sample, NbSupp et GeneID et grâce a un groupBy on construit la matrice
createMatrixComptage <- function(tableauDataCirc, JunctionCutOff){
  withProgress(message = "Create count matrix", value = 0.1, {
    tryCatch({
            incProgress(0.2, detail = "Get counts")
            data.shiny$TotalNumberOfGenesBeforeFiltering = length(unique(tableauDataCirc$GeneId)) - 1 # pour le charactère vide
            # On récupère les jonctions et graĉes au dcast on regroupe en une ligne celle qui ont même start, même end...
            JunctionData = dcast(tableauDataCirc, formula = Chromosome + Start + End + Strand + GeneId ~ Sample, value.var = "NbSuppReads", fill = 0)
            # On crée un coldta Spécifique à cette matrice
            JunctionColdata = data.shiny$coldata
            # JunctionColdata = subset.data.frame(tableauDataCirc, select = c("Sample", "Condition")) %>% distinct(Sample, Condition) %>% as.data.frame()
            JunctionData = JunctionData[,c("Chromosome", "Start", "End", "Strand", "GeneId", JunctionColdata$Sample)]
            ddsDE = NormalisationCirc(JunctionData[,c(JunctionColdata$Sample)], JunctionColdata)
            # On détermine les comptages normalisé et pour chaque jonction on regarde si elle vue à hauteur du CutOff au moins dans une condition mais on ne le fait pas direct sur les données car on doit anvoyé les comtages brut
            incProgress(0.4, detail = "Filter by junction cutoff")
            dataToTreat = counts(ddsDE, normalized=TRUE)
            CdnList = unique(JunctionColdata$Condition)
            Bool = apply(dataToTreat, 1, function(line){
            any(
                sapply(CdnList, function(cdn){
                all(line[which(JunctionColdata$Condition == cdn)] >= JunctionCutOff)
                })
            )
            })
            incProgress(0.7, detail = "Sum by genes")
            # Graĉe au Boolean on filtre les données. Attention: Le filtre ne modifie pas dataCirc, dans le summary, on les observe toujours
            FilteredJunctionData = JunctionData[Bool,]
            data.shiny$FilteredJunctionsFromCirc = FilteredJunctionData
            dataSub <- melt(FilteredJunctionData, id.vars = "GeneId", measure.vars = JunctionColdata$Sample, variable.name = "Sample", value.name = "NbSuppReads")
            cc <- dataSub %>% group_by(Sample, GeneId) %>% mutate(count = sum(NbSuppReads)) %>% as.data.frame() %>% subset.data.frame(select = c("Sample", "GeneId", "count")) %>% distinct(Sample, GeneId, count) %>% dcast(GeneId ~ Sample, fill = 0)
            cc <- cc[!is.na(cc[, 1]), ]
            cc <- data.frame(cc[, -1], row.names = cc[, 1])
            incProgress(1, detail = "Matrix done")
            return(cc)
    }, error = function(e){
        showNotification("Error in createMatrixComptage, can't create matrix from data !", type = "error")
        data.shiny$FilteredJunctionsFromCirc <- NULL
        return(NULL)
    })
})
  
}

# Getter
getCountMatrix <- reactive({
  req(input$mainEffectCirc)
  if (is.null(data.shiny$countMatrix)){
    data.shiny$countMatrix <- createMatrixComptage(data.shiny$dataCircRNAPrediction, input$JunctionCutOff)
    return(data.shiny$countMatrix)
  }else{
    return(data.shiny$countMatrix)
  }
})

#Créer la matrice de comptage
createCountMatrix <- observeEvent(input$reloadCountMatrix, {
  data.shiny$countMatrix <- createMatrixComptage(data.shiny$dataCircRNAPrediction, input$JunctionCutOff)
})


# Pour récupérer dans l'UI qu'elle est la condition 1
output$Cdn1 <- renderText({
  paste("Condition 1 : ", input$Condition_1_circ)
})
# Pour récupérer dans l'UI qu'elle est la condition 2
output$Cdn2 <- renderText({
  paste("Condition 2 : ", input$Condition_2_circ)
})


output$Number_of_genes_Raw_counts <- renderInfoBox({
    infoBox(
        "Reproductibility filter applied", paste(nrow(data.shiny$countMatrix), "genes"), icon = icon("filter"), 
        color = "purple"
    )
})


output$Number_of_genes_Filt_counts <- renderInfoBox({
    infoBox(
        "Basemean filter applied", paste(if (length(data.shiny$circRNANormalisedCounts) > 0) nrow(data.shiny$circRNANormalisedCounts) else nrow(data.shiny$countMatrix), "genes"), icon = icon("filter"), 
        color = "purple"
    )
})

output$Total_number_of_genes <- renderInfoBox({
    infoBox(
        "Total", paste(data.shiny$TotalNumberOfGenesBeforeFiltering, "genes"), icon = icon("database"),
        color = "purple"
    )
})


output$circDensplot <- renderHighchart({
  if (!is.null(data.shiny$countMatrix)){
    df = melt(as.matrix(data.shiny$countMatrix))
    colnames(df) = c("Gene", "Sample", "value")
    df$value = log(df$value + 1)
    hcboxplot(x = df$value, var = df$Sample) %>% hc_chart(type = "column") %>%
        hc_title(text = "Log(Count + 1) Distribution on circular RNA counts") %>% 
        hc_xAxis(title = "Samples") %>%
        hc_yAxis(title="Log(count +1)") %>% 
        hc_add_theme(hc_theme_smpl()) %>%
        hc_exporting(enabled = TRUE, filename = "CircDensityPlotCount") %>%
        hc_chart(zoomType = "xy")
  }
})


# Permet de représenté le graph des quantiles
# Pas besoin de action button pour le grpah vue qu'il y a deja un boutton pour la matrice
output$QuantilePlot <- renderHighchart({
  if (!is.null(data.shiny$countMatrix)){
    quantile =  apply(data.shiny$countMatrix, MARGIN = 2, function(x){return(quantile(x, seq(0,1,0.01)))})
    mdf = melt(quantile)
    colnames(mdf) = c("Quantile", "Condition", "Value")
    # On plot la courbe pour déterminer le coude
    hchart(mdf, "line", hcaes(x = Quantile, y = Value, group = Condition)) %>%  
        hc_title(text = "Quantile Plot") %>% 
        hc_yAxis(title="Quantile associated value") %>% 
        hc_xAxis(title = "Quantile") %>%
        hc_add_theme(hc_theme_smpl()) %>%
        hc_exporting(enabled = TRUE, filename = "Quantile Plot") %>%
        hc_chart(zoomType = "xy")
  }
  
})


# Permet de créer la table des Up et des Down
computeTableForUpDown <- function(DeSeqRes, upDownCutOff, pvAdjCutOff, isJunction){
  
  if (!is.null(DeSeqRes)){
        withProgress(message = "Get Up & Down", value = 0.1, {
            tryCatch({
                  incProgress(0.2, detail = "Applying cutoff")
                  RegulatedGeneUp = DeSeqRes[DeSeqRes$log2FoldChange > upDownCutOff & DeSeqRes$padj <= pvAdjCutOff, ]
                  RegulatedGeneDown = DeSeqRes[DeSeqRes$log2FoldChange < -upDownCutOff & DeSeqRes$padj <= pvAdjCutOff, ]
                  
                  BelowFoldChangeGeneUp = DeSeqRes[DeSeqRes$log2FoldChange > 0 & DeSeqRes$log2FoldChange < upDownCutOff & DeSeqRes$padj <= pvAdjCutOff, ]
                  BelowFoldChangeGeneDown = DeSeqRes[DeSeqRes$log2FoldChange < 0 & DeSeqRes$log2FoldChange > -upDownCutOff & DeSeqRes$padj <= pvAdjCutOff, ]
                  
                  BelowPvalueGeneUp = DeSeqRes[DeSeqRes$log2FoldChange > upDownCutOff & DeSeqRes$padj > pvAdjCutOff, ]
                  BelowPvalueGeneDown = DeSeqRes[DeSeqRes$log2FoldChange < -upDownCutOff & DeSeqRes$padj > pvAdjCutOff, ]
                  
                  BelowFC_PvalueGeneUp = DeSeqRes[DeSeqRes$log2FoldChange > 0 & DeSeqRes$log2FoldChange < upDownCutOff & DeSeqRes$padj > pvAdjCutOff, ]
                  BelowFC_PvalueGeneDown = DeSeqRes[DeSeqRes$log2FoldChange < 0 & DeSeqRes$log2FoldChange > -upDownCutOff & DeSeqRes$padj > pvAdjCutOff, ]
                  incProgress(0.6, detail = "Build tables")
                  tableau = data.frame(
                    Up = c(
                      nrow(RegulatedGeneUp), nrow(BelowFoldChangeGeneUp),nrow(BelowPvalueGeneUp),nrow(BelowFC_PvalueGeneUp),
                      sum(nrow(RegulatedGeneUp), nrow(BelowFoldChangeGeneUp), nrow(BelowPvalueGeneUp), nrow(BelowFC_PvalueGeneUp))
                    ),
                    Down = c(
                      nrow(RegulatedGeneDown), nrow(BelowFoldChangeGeneDown), nrow(BelowPvalueGeneDown), nrow(BelowFC_PvalueGeneDown),
                      sum(nrow(RegulatedGeneDown), nrow(BelowFoldChangeGeneDown),  nrow(BelowPvalueGeneDown), nrow(BelowFC_PvalueGeneDown))				
                    ),
                    Total = c(
                      sum(nrow(RegulatedGeneUp), nrow(RegulatedGeneDown)), sum(nrow(BelowFoldChangeGeneUp), nrow(BelowFoldChangeGeneDown)), sum(nrow(BelowPvalueGeneUp), nrow(BelowPvalueGeneDown)), sum(nrow(BelowFC_PvalueGeneUp), nrow(BelowFC_PvalueGeneDown)),
                      nrow(DeSeqRes)				
                    ),
                    row.names = c("Regulated", "Below FC", "Below padj", "Below FC and padj", "Total")
                    
                  )
                  if (isJunction){
                        data.shiny$UpDownJunction <- tableau
                  }else{
                        data.shiny$UpDownCirc <- tableau
                  }
                  incProgress(0.1, detail = "Done")
              }, error = function(e){
              showNotification("Error in computeTableForUpDown, can't access to table !", type = "error")
              return(data.frame())
            })
        })   
  }
  if (isJunction){
    return(data.shiny$UpDownJunction)
  }else{
    return(data.shiny$UpDownCirc)
  }
}

output$TableUpDownForCirc <- DT::renderDataTable(
    as.data.frame(t(data.shiny$UpDownCirc)),
    style = "bootstrap",
    rownames = TRUE,
    escape = FALSE,
    extensions = 'Buttons',
    options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)

)


output$PlotUpDownForCirc <- renderHighchart({
    if (!is.null(data.shiny$DEResultsOnlyOnCirc)){
        cc <- computeTableForUpDown(data.shiny$DEResultsOnlyOnCirc, input$CutOffFC_Circ, input$CutPadj_Circ, FALSE)[1:4, ]
        highchart() %>%
        hc_title(text = "Regulation resume for circular RNAs group by genes") %>% 
        hc_xAxis(categories = rownames(cc)) %>%
        hc_plotOptions(
            column = list(
            stacking = "normal"
            )
        ) %>%
        hc_add_series(
            name = "Up",
            data = cc$Up,
            type = "column"
        )%>% 
        hc_add_series(
            name = "Down",
            data = cc$Down,
            type = "column"
        ) %>%
        hc_exporting(enabled = TRUE, filename = "CircUpDown")
   }
})




# Fonction permettant de compute le test de Kolmigorov-Smirnov sur les jonctions associés aux n gènes 
computeKSTest <- function(){
  withProgress(message = "Compute KS test", value = 0.1, {
    tryCatch({
        incProgress(0.2, detail = "Get data")
        data <- data.shiny$dataForKSTest
        
        ProcessInfo <- data.shiny$dataProcessingInfo
        GeneList <- unique(data$GeneId)
        Conditions <- c(input$Condition_1_circ, input$Condition_2_circ)
        
        colnames(data)[-1] <- ProcessInfo[match(colnames(data)[-1], ProcessInfo$Sample, nomatch = 0), ]$Condition
        incProgress(0.5, detail = "Compute tests")
        Genes <- data$GeneId
        data <- sapply(Conditions, FUN = function(cdn){
        return(rowMeans(data[,c(which(colnames(data) == cdn))]))
        }) %>% as.data.frame()
        data$GeneId <- Genes
        test <- sapply(GeneList, FUN = function(x){
        KSData <- data[data$GeneId == x,]
        return(suppressWarnings(ks.test(x = KSData[,1], y = KSData[,2])$p.value))
        })
        res <- data.frame(Gene = GeneList, KS_p.value = test, KS_p.value_adj = p.adjust(test, "fdr"))
        incProgress(1, detail = "Done")
    }, error = function(e){
        showNotification("Error when computing KS Test !", type = "error")
        return(NULL)
    })
  })
  return(res)
}

launchKS_Test <- observeEvent(input$KSTestButt, {
  if (!is.null(data.shiny$dataForKSTest)){
    data.shiny$KSTestResults <- computeKSTest()
  }else{
    showNotification("To compute Kolmogorov Test, please first compute analysis on junctions to obtains normalised counts", type = "warning")
  }
})



output$KS_test_Table <- DT::renderDataTable(
    as.data.frame(data.shiny$KSTestResults),
    style = "bootstrap",
    rownames = TRUE,
    escape = FALSE,
    extensions = 'Buttons',
    options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)


output$Rep1_Circ <- renderUI({
  selectInput(inputId = "Rep1_Circ", label = "Replicat : 1", choices = unique(data.shiny$dataProcessingInfo$Replicat), selected = unique(data.shiny$dataProcessingInfo$Replicat)[1])
})

output$Rep2_Circ <- renderUI({
  selectInput(inputId = "Rep2_Circ", label = "Replicat : 2", choices = unique(data.shiny$dataProcessingInfo$Replicat), selected = unique(data.shiny$dataProcessingInfo$Replicat)[2])
})

EstimateDispersionperCondition <- function(){
  if (length(data.shiny$circRNANormalisedCounts) > 0 && !is.null(data.shiny$coldata) && !is.null(data.shiny$dataProcessingInfo)){
    withProgress(message = "Estimate dispersion", value = 0.1, {
        tryCatch({
            incProgress(0.2, detail = "Get data")
            coldata <- data.shiny$coldata
            ProcessInfo <- data.shiny$dataProcessingInfo
            incProgress(0.3, detail = "Build table by condition")
            res <- lapply(c(input$Condition_1_circ, input$Condition_2_circ), FUN = function(x){
            subTab = as.data.frame(log2(data.shiny$circRNANormalisedCounts[,as.vector(coldata[coldata$Condition == x,]$Sample)] + 1))
            colnames(subTab) = ProcessInfo[match(colnames(subTab), ProcessInfo$Sample, nomatch = 0), ]$Replicat
            subTab$Condition =  rep(x,nrow(subTab))
            subTab$GeneId = rownames(subTab)
            rownames(subTab) = c()
            return(subTab)
            })
            incProgress(0.7, detail = "Global table")
            table <- do.call("rbind", res)
            incProgress(1, detail = "Done")
        }, error = function(e){
            showNotification("Error in EstimateDispersionperCondition !", type = "error")
            return(NULL)
        })
    })
    return(table)
  }
}


output$DispersionPlot <- renderHighchart({
  data <- EstimateDispersionperCondition()
  if (!is.null(data)){
    data <- data[,c(input$Rep1_Circ, input$Rep2_Circ, "Condition", "GeneId")]
    colnames(data) = c("Rep1", "Rep2", "Condition", "GeneId")
    
    output$Rsquared_Cdn1_Circ <- renderText({
      cdn = input$Condition_1_circ
      dataForModel = data[data$Condition == cdn, ]
      model = lm(dataForModel, formula = Rep1 ~ Rep2)
      cc = summary(model)
      paste(cdn, " : ", cc$adj.r.squared)
    })
    
    output$Rsquared_Cdn2_Circ <- renderText({
      cdn = input$Condition_2_circ
      dataForModel = data[data$Condition == cdn, ]
      model = lm(dataForModel, formula = Rep1 ~ Rep2)
      cc = summary(model)
      paste(cdn, " : ", cc$adj.r.squared)
    })
    
    x <- c("GeneId : ")
    y <- c( sprintf("{point.%s}", c("GeneId")))
    tltip <- tooltip_table(x, y)
    hc <- highchart()
    hchart(data, "scatter", hcaes(x = Rep1, y = Rep2, group = Condition))  %>%
      hc_tooltip(useHTML = TRUE, headerFormat = "", pointFormat = tltip) %>% 
      hc_xAxis(title = paste("Log2(count +1)", input$Rep1_Circ)) %>%
      hc_yAxis(title = paste("Log2(count +1)", input$Rep2_Circ))
  }
})


# Permet de compute l'analyse uniquement sur les comptages des circulaires
computeDEAnalysisOnlyOnCirc <- function(circCountMatrix, coldata){
  if (!is.null(circCountMatrix) && !is.null(coldata)){
    withProgress(message = "DE Analysis on genes", value = 0.1, {
        tryCatch({
            incProgress(0.2, detail = "Get counts")
            circCountMatrix <- circCountMatrix[, c(coldata$Sample)]

            incProgress(0.4, detail = "Applying baseMean filter")
            # First pour déterminer les gènes à exclure avec la normalisation choisie
            dds <- NormalisationCirc(circCountMatrix, coldata)
            filteredCountMatrix <- circCountMatrix[rowMeans(counts(dds, normalized=TRUE)) >= input$BasemeanCutCirc, ]
            
            incProgress(0.6, detail = "Normalise counts")
            ddsDE = NormalisationCirc(filteredCountMatrix, coldata)

            incProgress(0.7, detail = "Launch DE Analysis")
            if(input$DE_Test_circ == 2){
                batchDesign = if(!is.null(input$batchEffectCirc)) paste("~", paste(input$batchEffectCirc, collapse = " + ")) else "~ 1"
                print(batchDesign)
                ddsDE <- DESeq(ddsDE, test = "LRT", reduced = formula(batchDesign))
            }else{
                ddsDE <- DESeq(ddsDE, test = "Wald")
                print(design(ddsDE))
            }
            incProgress(0.8, detail = "Extract normalised counts")
            data.shiny$circRNANormalisedCounts <- counts(ddsDE, normalized=TRUE)
            
            incProgress(0.9, detail = "Parse results")
            data.shiny$DEResultsOnlyOnCirc <- results(ddsDE, pAdjustMethod = "fdr", contrast = c(input$mainEffectCirc, input$Condition_2_circ, input$Condition_1_circ))
            output$normFactorCircByGenes <- renderTable(data.frame(Sample = names(sizeFactors(ddsDE)), Factor = sizeFactors(ddsDE)), label = "Normalisation Factor")
            
            incProgress(1, detail = "Done")
            }, error = function(e){
            
            showNotification("Error in computeDEAnalysisOnlyOnCirc !", type = "error")
            data.shiny$DEResultsOnlyOnCirc <- NULL       
            })
    })
    withProgress(message = "Treat Na values", value = 0.1, {
        incProgress(0.2, detail = "On log2FC")
        if (sum(is.na(data.shiny$DEResultsOnlyOnCirc$log2FoldChange)) > 0){
            data.shiny$DEResultsOnlyOnCirc[is.na(data.shiny$DEResultsOnlyOnCirc$log2FoldChange), ]$log2FoldChange  <- 0
        }
        incProgress(0.5, detail = "On p-values")
        if (sum(is.na(data.shiny$DEResultsOnlyOnCirc$padj)) > 0){
            data.shiny$DEResultsOnlyOnCirc[is.na(data.shiny$DEResultsOnlyOnCirc$padj), ]$padj  <- 1
        }
    incProgress(1, detail = "Done")    
    })
  }
}

launchDEAnalysisOnlyOnCirc <- observeEvent(input$DEAnalysisOnlyCircButt, {
  computeDEAnalysisOnlyOnCirc(getCountMatrix(), data.shiny$coldata)
  getGeneSymbolForCirc()
  getGeneBiotypeForCirc()
  
})

getGeneSymbolForCirc <- function(){
  if (!is.null(data.shiny$DEResultsOnlyOnCirc)){
    sub = left_join(x=(data.frame(GeneId = rownames(data.shiny$DEResultsOnlyOnCirc))), y=unique((subset.data.frame(data.shiny$dataCircRNAPrediction, select = c("GeneId", "GeneSymbol")))), by =  'GeneId', copy=FALSE)
    data.shiny$DEResultsOnlyOnCirc$gene_symbol = sub$GeneSymbol
  }
}

getGeneBiotypeForCirc <- function(){
  if (!is.null(data.shiny$DEResultsOnlyOnCirc)){
    sub = left_join(x=(data.frame(GeneId = rownames(data.shiny$DEResultsOnlyOnCirc))), y=unique((subset.data.frame(data.shiny$dataCircRNAPrediction, select = c("GeneId", "GeneBioType")))), by =  'GeneId', copy=FALSE)
    data.shiny$DEResultsOnlyOnCirc$GeneBioType = sub$GeneBioType
  }
}


# Permet de créer le graph de densité des comptage pour les comptage en ARN circulaire sur les gènes
# Pas besoin de action button pour le grpah vue qu'il y a deja un boutton pour la matrice
output$circDensplotNormalised <- renderHighchart({
  if (!is.null(data.shiny$circRNANormalisedCounts)){
    df = melt(as.matrix(data.shiny$circRNANormalisedCounts))
    colnames(df) = c("Gene", "Sample", "value")
    df$value = log(df$value + 1)
    hcboxplot(x = df$value, var = df$Sample) %>% hc_chart(type = "column") %>%
      hc_title(text = "Log(Count + 1) Distribution on circular RNA counts for normalised filtered circular RNA counts") %>% 
      hc_xAxis(title = "Samples") %>%
      hc_yAxis(title="Log(count +1)") %>% 
      hc_add_theme(hc_theme_smpl()) %>%
      hc_exporting(enabled = TRUE, filename = "CircDensityPlotCount") %>%
      hc_chart(zoomType = "xy")
  }
})

output$DEResultsTableOnlyOnCirc <- DT::renderDataTable(
  as.data.frame(data.shiny$DEResultsOnlyOnCirc),
  style="bootstrap",
  rownames = TRUE,
  escape = FALSE,
  extensions = 'Buttons',
  options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)

computeDEAnalysisOnCircularJunctions = function(CountsOnJunctions, coldata){
  if (!is.null(CountsOnJunctions) && !is.null(coldata)){
    withProgress(message = "DE Analysis on junctions", value = 0.1, {
        tryCatch({
                incProgress(0.2, detail = "Get counts")
                JunctionsCircCountMatrix = CountsOnJunctions[, c(coldata$Sample)]
                rownames(JunctionsCircCountMatrix) = CountsOnJunctions %>% rowwise() %>% transmute(Id = paste(Chromosome, Start, End, Strand, GeneId, collapse = "", sep = "_")) %>% unlist()
                
                incProgress(0.4, detail = "Applying baseMean filter")
                dds <- NormalisationCirc(JunctionsCircCountMatrix, coldata)
                # First pour déterminer les gènes à exclure avec la normalisation de DESeq2
                bool <- rowMeans(counts(dds, normalized=TRUE)) >= input$BasemeanCutJunctions
                filteredCountMatrix <- JunctionsCircCountMatrix[bool, ]
                
                incProgress(0.6, detail = "Normalise counts")
                # Second pour réaliser la DE
                ddsDE <- NormalisationCirc(filteredCountMatrix, coldata)
                
                incProgress(0.7, detail = "Launch DE Analysis")
                if (input$DE_Test_circ == 2){
                    batchDesign <- if (!is.null(input$batchEffectCirc)) paste("~", paste(input$batchEffectCirc, collapse = " + ")) else "~ 1"
                    print(batchDesign)
                    ddsDE <- DESeq(ddsDE, test = "LRT", reduced = formula(batchDesign))
                }else{
                    ddsDE <- DESeq(ddsDE, test = "Wald")
                    print(design(ddsDE))
                }
                
                incProgress(0.8, detail = "Extract normalised counts")
                output$normFactorCircJunctions <- renderTable(data.frame(Sample = names(sizeFactors(ddsDE)), Factor = sizeFactors(ddsDE)), label = "Normalisation Factor")
                data.shiny$JunctionNormalisedCounts = counts(ddsDE, normalized=TRUE)
                
                data.shiny$FilteredJunctionsFromCircDescriptions = subset.data.frame(CountsOnJunctions[bool, ], select = c("Chromosome", "Start", "End", "Strand", "GeneId"))
                # On récupère les lignes où le gène et au moins associé à 2 jonctions pour effectuer le test de Fisher il nous faut au moins deux jonctions
                data.shiny$dataForKSTest = cbind(select_(CountsOnJunctions[bool, ], "GeneId"), data.shiny$JunctionNormalisedCounts)  %>% group_by(GeneId) %>% filter(n() >=2 & GeneId != "")
                incProgress(0.9, detail = "Parse results")
                data.shiny$DEResultsOnlyCircularJunctions = results(ddsDE, pAdjustMethod = "fdr", contrast = c(input$mainEffectCirc, input$Condition_2_circ, input$Condition_1_circ))
        
                incProgress(1, detail = "Done")
        }, error = function(e){
            
            showNotification("Error in computeDEAnalysisOnCircularJunctions !", type = "error")
            data.shiny$DEResultsOnlyCircularJunctions <- NULL
            
        })
    incProgress(1, detail = "Done")
    })
    withProgress(message = "Treat Na values", value = 0.1, {
        incProgress(0.2, detail = "On log2FC")
        if(sum(is.na(data.shiny$DEResultsOnlyCircularJunctions$log2FoldChange)) > 0){
            data.shiny$DEResultsOnlyCircularJunctions[is.na(data.shiny$DEResultsOnlyCircularJunctions$log2FoldChange), ]$log2FoldChange  <- 0
        }
        incProgress(0.5, detail = "On p-values")
        if(sum(is.na(data.shiny$DEResultsOnlyCircularJunctions$padj)) > 0){
            data.shiny$DEResultsOnlyCircularJunctions[is.na(data.shiny$DEResultsOnlyCircularJunctions$padj), ]$padj  <- 1
        }
        incProgress(1, detail = "Done")
    })
  }
}


lauchDEOnJunctionsCirculars <- observeEvent(input$DEAnalysisOnlyJunctionsButt, {
  if(length(data.shiny$FilteredJunctionsFromCirc) > 0){
    computeDEAnalysisOnCircularJunctions(data.shiny$FilteredJunctionsFromCirc, data.shiny$coldata)
  }
})

output$DEResultsTableOnlyOnJunctionsCirc <- DT::renderDataTable(
    as.data.frame(data.shiny$DEResultsOnlyCircularJunctions),
    style="bootstrap",
    rownames = TRUE,
    escape = FALSE,
    extensions = 'Buttons',
    options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)

output$TableUpDownForJunctionsCirc <- DT::renderDataTable(

    as.data.frame(t(data.shiny$UpDownJunction)),
    style="bootstrap",
    rownames = TRUE,
    escape = FALSE,
    extensions = 'Buttons',
    options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)


output$PlotUpDownForJunctionsCirc <- renderHighchart({
    if(!is.null(data.shiny$DEResultsOnlyCircularJunctions)){
        cc = computeTableForUpDown(data.shiny$DEResultsOnlyCircularJunctions, input$CutOffFC_Circ, input$CutPadj_Circ, TRUE)[1:4, ]
        highchart() %>%
        hc_title(text = "Regulation resume for circular RNAs junctions") %>%
        hc_xAxis(categories = rownames(cc)) %>%
        hc_plotOptions(
            column = list(
            stacking = "normal"
            )
        ) %>%
        hc_add_series(
            name = "Up",
            data = cc$Up,
            type = "column"
        )%>%
        hc_add_series(
            name = "Down",
            data = cc$Down,
            type = "column"
        ) %>%
        hc_exporting(enabled = TRUE, filename = "JunctionsCircUpDown")
  }
})

output$CountsNormalisedForCircularOnGenes <- DT::renderDataTable(

    data.frame(as.matrix(GeneId = rownames(data.shiny$circRNANormalisedCounts), data.shiny$circRNANormalisedCounts, rownames = FALSE)),
    style="bootstrap",
    rownames = TRUE,
    escape = FALSE,
    filter = "top",
    extensions = 'Buttons',
    options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)

output$RawCountsForJunctionsThatBuildGenes <- DT::renderDataTable(
 
    as.data.frame(data.shiny$FilteredJunctionsFromCirc),
    style="bootstrap",
    rownames = TRUE,
    escape = FALSE,
    filter = "top",
    extensions = 'Buttons',
    options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)

)

output$CountsNormalisedForCircularJunctions <- DT::renderDataTable(

    data.frame(data.shiny$FilteredJunctionsFromCircDescriptions, data.shiny$JunctionNormalisedCounts),
    style="bootstrap",
    rownames = TRUE,
    escape = FALSE,
    filter = "top",
    extensions = 'Buttons',
    options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)

)
