library(patchwork)
library(shiny)
library(Matrix)
library(shinydashboard)
library(ggplot2)
library(shinyWidgets)
library(shinyjs)
library(shinyFiles)
library(ggrepel)
library(stringr)
library(dplyr)
library(Seurat)
library(hdf5r)
library(readr)
library(clustree)
library(tidyverse)
library(fontawesome)
library(cowplot)
library(Signac)
library(DT)
library(data.table)
library(plotly)
library(colourpicker)
library(httr)
library(jsonlite)
library(monocle3)
library(SingleCellExperiment)
library(SeuratWrappers)


# Importer les données
options(shiny.maxRequestSize = 2000*1024^2)  # 2000 MB


# Source ui files:


# Source server files:
source("single_dataset_server.R")
source("multiple_datasets_server.R")
source("trajectory_server.R")

source("main_ui.R")



# Main Server Function
server <- function(input, output, session) {
  single_dataset_server(input, output, session)
  multiple_datasets_server(input, output, session)
  trajectory_server(input, output, session)
}

# Exécuter l'application
options(shiny.launch.browser = TRUE)

shinyApp(ui = ui, server = server)
