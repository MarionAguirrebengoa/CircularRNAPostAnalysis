
# Permet de créer le graph de densité des comptage globaux sur les gènes
# Pas beosin  de action button vue que sa sera jamais filtrer
output$geneDensplot <- renderHighchart({
  if (!is.null(data.shiny$dataGeneCount)){
    df <- melt(as.matrix(data.shiny$dataGeneCount))
    colnames(df) <- c("Gene", "Sample", "value")
    df$value <- log(df$value + 1)
        hcboxplot(x = df$value, var = df$Sample) %>% hc_chart(type = "column") %>%
        hc_title(text = "Log(Count + 1) Distribution on Gene counts") %>% 
        hc_xAxis(title = "Samples") %>%
        hc_yAxis(title="Log(count +1)") %>% 
        hc_add_theme(hc_theme_smpl()) %>%
        hc_exporting(enabled = TRUE, filename = "GeneDensityPlotCount") %>%
        hc_chart(zoomType = "xy")
  }
})


output$mainEffect = renderUI({
  selectInput(inputId = "mainEffect", label = "main effect to test", choices = colnames(data.shiny$coldata), selectize = TRUE, multiple = FALSE)
})

output$batchEffect = renderUI({
  selectInput(inputId = "batchEffect", label = "batch effect", choices = colnames(data.shiny$coldata), selectize = TRUE, multiple = TRUE)
})

output$Condition_1_gene = renderUI({
if (!is.null(input$mainEffect)){
      selectInput(inputId = "Condition_1_gene", label = "Condition 1", choices = data.shiny$coldata[[input$mainEffect]], selectize = TRUE, multiple = FALSE)
}
})

output$Condition_2_gene = renderUI({
  if (!is.null(input$mainEffect)){
      selectInput(inputId = "Condition_2_gene", label = "Condition 2", choices = data.shiny$coldata[[input$mainEffect]], selectize = TRUE, multiple = FALSE)
  }
})

NormalisationGene = function(CountMatrix, coldata){
  dds <- NULL
  withProgress(message = "Normalise counts", value = 0.1, {
        tryCatch({
            incProgress(0.2, detail = "Get design to prepare object")
            Design = paste("~", paste( paste(input$mainEffect, collapse = " + ") , paste(input$batchEffect, collapse = " + "), sep = ifelse(!is.null(input$batchEffect), yes = " + ", no = "")))
            # HouseKeeping genes
            if (input$NormMethod == 1){
                incProgress(0.4, detail = "House Keeping genes")
                # On prépare l'objet DeSeq
                dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = coldata, design = formula(Design))
                # On fait la normalisation avec ls HKG
                dds <- estimateSizeFactors(object = dds, controlGenes=(rownames(CountMatrix) %in% input$HouseKeepingGenes))
            }
            # CPM normalisation
            else if (input$NormMethod == 2){
                incProgress(0.4, detail = "CPM")
                NbMappedReads <- colSums(CountMatrix)
                CPMFactor <- NbMappedReads/input$normFactCPM
                # On prépare l'objet DeSeq
                dds = DESeqDataSetFromMatrix(countData = CountMatrix, colData = coldata, design = formula(Design))
                # On met le sizeFactor en CPM
                sizeFactors(dds) = CPMFactor
            }else if (input$NormMethod == 3){
                incProgress(0.4, detail = "Deseq")
                # On fait un objet classique
                dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = coldata, design = formula(Design))
                dds <- estimateSizeFactors(dds)
            }else{
                incProgress(0.4, detail = "RPM")
                NbMappedReads <- data.shiny$dataProcessingInfo[match(coldata$Sample, data.shiny$dataProcessingInfo$Sample),]$NbReadsMappedBySTAR
                CPMFactor <- NbMappedReads/input$normFactCPM
                dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = coldata, design = formula(Design))
                sizeFactors(dds) <- CPMFactor
            }
            incProgress(1, detail = "Normalisation done")
        }, error = function(e){
            showNotification("Error in NormalisationGene, can't normalise data !", type = "error")
            showNotification(e)
            dds <- NULL
        })
  })
  return(dds)
 
}

output$HouseKeepingGenes <- renderUI({
  if (!is.null(data.shiny$dataGeneCount)){
    withProgress(message = "get HKG", value = 0.1, {
        GeneIdList <- rownames(data.shiny$dataGeneCount)
        incProgress(1, detail = "Done")
    })
    selectInput(inputId = "HouseKeepingGenes", label = "House Keeping Genes", choices = GeneIdList, selectize = TRUE, multiple = TRUE)
  }
})



output$Number_of_genes_Filt_counts_Gene <- renderInfoBox({
  if (!is.null(data.shiny$geneNormalisedCounts)){
    infoBox(
        "BaseMean filter applied", paste(if(length(data.shiny$geneNormalisedCounts) > 0) nrow(data.shiny$geneNormalisedCounts) else data.shiny$TotalNumberOfGenesFromHTSeqCounts, "genes"), icon = icon("filter"),
        color = "purple"
    )
  }
})

output$Total_number_of_genes_HTSeq <- renderInfoBox({
  if (!is.null(data.shiny$TotalNumberOfGenesFromHTSeqCounts)){
    infoBox(
        "Total", paste(data.shiny$TotalNumberOfGenesFromHTSeqCounts, "genes"), icon = icon("database"), 
        color = "purple"
    )
  }
})

# Permet de compute l'analyse uniquement sur les comptages des gènes
computeDEAnalysisOnlyOnGene <- function(geneCountMatrix, coldata){
  if (!is.null(geneCountMatrix) && !is.null(coldata)){
    withProgress(message = "DE Analysis on genes", value = 0.1, {
        tryCatch({
            incProgress(0.2, detail = "Get counts")
            incProgress(0.4, detail = "Applying baseMean filter")
            # On normalise une première fois pour le filtre, ici pour le filtre le design n'est pas important puisqu'il ne sera pas utilisé mais on évite de dupliquer la fonction
            dds = NormalisationGene(geneCountMatrix, coldata)
            filteredCountMatrix = geneCountMatrix[rowMeans(counts(dds, normalized=TRUE)) >= input$BasemeanCutGene, ]
            
            incProgress(0.6, detail = "Normalise counts")
            # Puis on lance l'analyse
            ddsDE = NormalisationGene(filteredCountMatrix, coldata)
            
            incProgress(0.7, detail = "Launch DE Analysis")
            # On détermine le test à appliquer selon si on est en simple ou multiple factor
            if(input$DE_Test_gene == 2){
                batchDesign = if(!is.null(input$batchEffect)) paste("~", paste(input$batchEffect, collapse = " + ")) else "~ 1"
                print(batchDesign)
                ddsDE <- DESeq(ddsDE, test = "LRT", reduced = formula(batchDesign))
            }else{
                ddsDE <- DESeq(ddsDE, test = "Wald")
            }
            incProgress(0.8, detail = "Extract normalised counts")
            data.shiny$geneNormalisedCounts <- counts(ddsDE, normalized=TRUE)

            incProgress(0.9, detail = "Parse results")
            data.shiny$DEResultsOnlyOnGene <- results(ddsDE, pAdjustMethod = "fdr", contrast = c(input$mainEffect, input$Condition_2_gene, input$Condition_1_gene))
            data.shiny$normFactorForGenes <- sizeFactors(ddsDE)
            output$normFactorGenes <- renderTable(data.frame(Sample = names(data.shiny$normFactorForGenes), Factor = data.shiny$normFactorForGenes), label = "Normalisation Factor")
            
            incProgress(1, detail = "Done")
        }, error = function(e){
            showNotification("Error in computeDEAnalysisOnlyOnGene !", type = "error")
            data.shiny$DEResultsOnlyOnGene <- NULL
        })
    })
    withProgress(message = "Treat Na values", value = 0.1, {
        incProgress(0.2, detail = "On log2FC")
        if (sum(is.na(data.shiny$DEResultsOnlyOnGene$log2FoldChange)) > 0){
            data.shiny$DEResultsOnlyOnGene[is.na(data.shiny$DEResultsOnlyOnGene$log2FoldChange), ]$log2FoldChange  <- 0
        }
        incProgress(0.5, detail = "On p-values")
        if (sum(is.na(data.shiny$DEResultsOnlyOnGene$padj)) > 0){
            data.shiny$DEResultsOnlyOnGene[is.na(data.shiny$DEResultsOnlyOnGene$padj), ]$padj  <- 1
        }
        incProgress(1, detail = "Done")
    })
  }
}



# Permet de créer la table des Up et des Down
computeTableForUpDownGene = function(DeSeqRes, upDownCutOff, pvAdjCutOff){

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
                data.shiny$UpDownGene = data.frame(
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
                incProgress(0.1, detail = "Done")
        }, error = function(e){
          showNotification("Error in computeTableForUpDownGene, process table !", type = "error")
          return(data.frame())
        })
      })
  }

  return(data.shiny$UpDownGene)
}


output$TableUpDownForGene <- DT::renderDataTable(
  as.data.frame(t(data.shiny$UpDownGene)),
  style="bootstrap",
  rownames = TRUE,
  escape = FALSE,
  extensions = 'Buttons',
  options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)


output$PlotUpDownForGene <- renderHighchart({
if (!is.null(data.shiny$DEResultsOnlyOnGene)){
        cc <- computeTableForUpDownGene(data.shiny$DEResultsOnlyOnGene, input$CutOffFC_Gene, input$CutPadj_Gene)[1:4,]
        highchart() %>%
        hc_title(text = "Regulation resume for Genes") %>% 
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
        hc_exporting(enabled = TRUE, filename = "GeneUpDown")
 }
})



laucnDEAnalysisOnlyOnGene <- observeEvent(input$DEAnalysisOnlyGeneButt, {
  computeDEAnalysisOnlyOnGene(data.shiny$dataGeneCount, data.shiny$coldata)
  if (input$Organism != "None" && input$Idtype != "None"){
    getGeneSymbolForGene()
    getGeneBioTypeForGene()
  }
})


getGeneSymbolForGene <- function(){
  withProgress(message = "Get gene symbols", value = 0.1, {
    tryCatch({
        incProgress(0.2, detail = "Get genes ids")
        data.shiny$ensembl <- useDataset(dataset = input$Organism, mart = data.shiny$ensembl)
        myDataGeneID <- as.data.frame(sapply(strsplit(x = rownames(data.shiny$DEResultsOnlyOnGene), split = "[.]"), '[', 1));
        colnames(myDataGeneID) <- "geneID"
        #biomaRt
        incProgress(0.4, detail = "BiomaRt query")
        biomartGeneSymbol <- biomaRt::getBM(attributes=c(input$Idtype, 'hgnc_symbol'), filters = input$Idtype, values = myDataGeneID$geneID, mart = data.shiny$ensembl) 
        biomartGeneSymbol <- data.frame(lapply(biomartGeneSymbol, as.character), stringsAsFactors = TRUE);
        colnames(biomartGeneSymbol) <- c("geneID", "hgnc_symbol")
        incProgress(0.8, detail = "Parse information")
        #filterResults car biomaRt peut renvoyer plusieurs résultata pour un même id, on prend le premier
        biomartGeneSymbol <- biomartGeneSymbol[!duplicated(biomartGeneSymbol$geneID), ]
        output <- dplyr::left_join(myDataGeneID, as.data.frame(biomartGeneSymbol), by = "geneID")
        #On associe les identifiants
        data.shiny$DEResultsOnlyOnGene$gene_symbol = output$hgnc_symbol;
        incProgress(1, detail = "Done")
        }, error = function(e){
            showNotification("Error in getGeneSymbolForGene, Check if BiomaRt is down ?", type = "error")
            data.shiny$DEResultsOnlyOnGene$gene_symbol <- NULL
        })
  })
}


getGeneBioTypeForGene <- function(){
    withProgress(message = "Get gene biotypes", value = 0.1, {
        tryCatch({
            incProgress(0.2, detail = "Get genes ids")
            myDataGeneID <- as.data.frame(sapply(strsplit(x = rownames(data.shiny$DEResultsOnlyOnGene), split = "[.]"), '[', 1));
            colnames(myDataGeneID) <- "geneID"
            incProgress(0.4, detail = "BiomaRt query")
            biomartGeneBiotype <- biomaRt::getBM(attributes=c(input$Idtype, 'gene_biotype'), filters = input$Idtype, values = myDataGeneID$geneID, mart = data.shiny$ensembl) 
            biomartGeneBiotype <- data.frame(lapply(biomartGeneBiotype, as.character), stringsAsFactors = TRUE);
            colnames(biomartGeneBiotype) <- c("geneID", "gene_biotype");
            incProgress(0.8, detail = "Parse information")
            #filterResults car biomaRt peut renvoyer plusieurs résultata pour un même id, on prend le premier
            biomartGeneBiotype <- biomartGeneBiotype[!duplicated(biomartGeneBiotype$geneID), ]
            output <-  dplyr::left_join(myDataGeneID, as.data.frame(biomartGeneBiotype), by = "geneID")
            #On associe les identifiants
            data.shiny$DEResultsOnlyOnGene$gene_biotype = output$gene_biotype;
            incProgress(1, detail = "Done")
            }, error = function(e){
                showNotification("Error in getGeneBioTypeForGene, Check if BiomaRt is down ?", type = "error")
                data.shiny$DEResultsOnlyOnGene$gene_biotype <- NULL
            })
    })
}


# Permet de créer le graph de densité des comptage globaux sur les gènes
# Pas beosin  de action button vue que sa sera jamais filtrer
output$geneDensplotNormalised <- renderHighchart({
  if (!is.null(data.shiny$geneNormalisedCounts)){
    df = melt(as.matrix(data.shiny$geneNormalisedCounts))
    colnames(df) = c("Gene", "Sample", "value")
    df$value = log(df$value + 1)
    hcboxplot(x = df$value, var = df$Sample) %>% hc_chart(type = "column") %>%
      hc_title(text = "Log(Count + 1) Distribution on normalised Gene counts") %>% 
      hc_xAxis(title = "Samples") %>%
      hc_yAxis(title="Log(count +1)") %>% 
      hc_add_theme(hc_theme_smpl()) %>%
      hc_exporting(enabled = TRUE, filename = "GeneDensityPlotCount") %>%
      hc_chart(zoomType = "xy")
  }
})


output$DEResultsTableOnlyOnGene <- DT::renderDataTable(
  as.data.frame(data.shiny$DEResultsOnlyOnGene),
  style="bootstrap",
  rownames = TRUE,
  escape = FALSE,
  extensions = 'Buttons',
  options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)



output$Organism <- renderUI({
  data <- listDatasets(data.shiny$ensembl)
  selectInput("Organism", "Choose organism", choices = c("", setNames(data$dataset, data$description), "None"="None"), selected =  "None", multiple = FALSE)
})

output$GeneNormalisedCounts <- DT::renderDataTable(
  data.frame(GeneId = rownames(data.shiny$geneNormalisedCounts), as.matrix(data.shiny$geneNormalisedCounts, rownames = FALSE)),
  style="bootstrap",
  rownames = TRUE,
  escape = FALSE,
  filter = "top",
  extensions = 'Buttons',
  options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'), lengthMenu = list(c(25, 50, 100,  -1), c('25', '50', '100', 'All Candidates')), pageLength = 15)
)
