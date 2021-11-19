tabItem(tabName = "DataImport",
        #busyIndicator("Calculation In progress", wait = 0),
        fluidRow(
          headerPanel("Import RDS files"),
          br(),
          br(),
          h3("Please always refer Circular RNA Prediction files before Gene read counts"),
          fluidRow(
                  column(6,
                         box(title = "Circular RNA Prediction files", width = 12, status = "primary",
                               fileInput("circRNAPred", "Import Circular RNA Prediction file",
                                         multiple = FALSE,
                                         accept = c(".rds"))
                         )
                  ),
                  column(6,
                        box(title = "Gene read counts", width = 12, status = "primary",
                               fileInput("GeneReadCount", "Import gene read count file",
                                         multiple = FALSE,
                                         accept = c(".rds"))
                               )
                )
        ),
        fluidRow(
          column(6, 
                 box(title="Processing Info files", width = 12, status = "primary",
                       fileInput("ProcessInfo", "Import Processing Info file",
                                 multiple = FALSE,
                                 accept = c(".rds"))
                 )
          ),
          column(6, 
                 box(title="Phenotype Informations", width = 12, status = "primary",
                     fileInput("PhenotypeInfo", "Import phenotype information",
                               multiple = FALSE,
                               accept = c(".csv"))
                 )
          )
        )
        )
)