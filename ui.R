#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(shinyjs)
library(shinybusy)
#library(shinysky)

source("theme.R", local = TRUE)$value

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
  includeCSS("styles.css"),
  useShinyjs(),
    dashboardPage(
      dashboardHeader(title = "Circular RNA Post Analysis", titleWidth = 300
      ),
      dashboardSidebar(width = 300, sidebarMenu(
        menuItem(text = "Data Import", icon = icon("archive"), tabName = "DataImport"),
        menuItem(text = "Summary", icon = icon("info-circle "), tabName = "Summary"),
        menuItem(text = "Filter data", icon = icon("filter"), tabName = "Filter"),
        menuItem(text = "Primary Statistics", icon = icon("percent"), tabName = "PrimaryStat"),
        menuItem(text = "Differential Expression Analysis", icon = icon("signal"), tabName = "DEAnalysis", startExpanded = FALSE,
                 menuSubItem("DE Analysis on Circular", tabName = "DEOnCirc"),
                 menuSubItem("DE Analysis on Gene", tabName = "DEOnGene"),
                 menuSubItem("Common DE Analysis on Circular And Gene", tabName = "DEOnCircGene")
                 )
      )),

      dashboardBody(
        theme_purple_gradient,
        tabItems(
          source(file.path("ui", "load_dataUI.R"), local = TRUE)$value,
          source(file.path("ui", "SummaryUI.R"), local = TRUE)$value,
          source(file.path("ui", "filterUI.R"), local = TRUE)$value,
          source(file.path("ui", "PrimaryStatUI.R"), local = TRUE)$value,
          source(file.path("ui", "DEAnalysisUICirc.R"), local = TRUE)$value,
          source(file.path("ui", "DEAnalysisUIGene.R"), local = TRUE)$value,
          source(file.path("ui", "DEAnalysisUICircGene.R"), local = TRUE)$value
        )
      )
    )
  )
)