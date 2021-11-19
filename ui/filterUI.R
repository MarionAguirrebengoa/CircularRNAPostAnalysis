
tabItem(tabName = "Filter",
        # busyIndicator("Calculation In progress", wait = 0),
        headerPanel("Filter"),
        fluidRow(align = "center",
                column(4, offset = 4,
                box(status = "warning", width=12,
                    column(6,
                            actionButton("filterButton", "Filter")
                    ),
                    column(6,
                            actionButton("ResetButton", "Reset")
                    )
                )
                )
        ),
        br(),
        br(),
        br(),
        fluidRow(align = "center",
                 column(6,
                     box(title = "Cross - references", status = "info", width = 12,
                           uiOutput("Genes"),
                           uiOutput("GeneId"),
                           uiOutput("GeneBiotypes")
                     )
                 ),
                 column(6,
                        box(title = "Support", status = "info", width = 12,
                            uiOutput("JunctionsReads"),
                            uiOutput("SupportedMethods")
                        )
                 )
        ),
        fluidRow(align = "center",
                  column(6,
                         fluidRow(
                           box(title = "General", status = "info", width = 12,
                                 uiOutput("Replicate"),
                                 uiOutput("Conditions"),
                                 uiOutput("Samples")
                           )
                         ),
                        fluidRow(
                            box(title = "Repeat regions", status = "info", width = 12,
                                column(6,
                                    textOutput("filechosen")
                                ),
                                column(6,
                                    actionButton("filechoose", label = "Import .bed repeat region file")
                                ),
                                    actionButton("removeRepeatRegions", label = "Remove repeat regions")
                            )
                        )
                  ),
                  column(6,
                         fluidRow(align = "center",
                             box(title = "Structural", status = "info", width = 12,
                                 column(6,
                                     uiOutput("Chromosomes"),
                                     box(title = "Length", status = "primary", width = 12,
                                           uiOutput("LengthMin"),
                                           uiOutput("LengthMax")
                                     )
                                 ),
                                 column(6,
                                 uiOutput("Signal"),
                                 uiOutput("Classes"),
                                 uiOutput("ExactMatch"),
                                 uiOutput("isChimera")
                                 ),
                                 box(title = "Merging information", status = "primary", width = 12,
                                    column(6,
                                    checkboxInput("globalMerged", "use global merged info (no matter sample)", value = FALSE, width = NULL),
                                    uiOutput("useMerge")
                                    ),
                                    column(6,
                                    checkboxInput("globalUnmerged", "use global unmerged info (no matter sample)", value = FALSE, width = NULL),
                                    uiOutput("useUnmerge")
                                    )
                                 )
                             )
                         )
                  )
        )
        
        
        
)
