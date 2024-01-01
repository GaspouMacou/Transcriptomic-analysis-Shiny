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
library(SeuratDisk)
library(data.table)
library(plotly)
library(colourpicker)
library(httr)
library(reticulate)
library(anndata)
library(SeuratObject)
#library(SeuratData)
library(jsonlite)
#library(monocle3)
#library(SeuratWrappers)
#Genomes
#library(EnsDb.Hsapiens.v86)
#library(biovizBase)



# Importer les données
options(shiny.maxRequestSize = 20000*1024^2)  # 2000 MB


# Source ui files:


# Source server files:
source("single_dataset_server.R")
source("multiple_datasets_server.R")
source("trajectory_server.R")
source("atac_server.R")
source("converter_server.R")
source("main_ui.R")


server <- function(input, output, session) {
  single_dataset_server(input, output, session)
  multiple_datasets_server(input, output, session)
  atac_server(input, output, session)
  trajectory_server(input, output, session)
  converter_server(input, output, session)

}

# Exécuter l'application
options(shiny.launch.browser = TRUE)

shinyApp(ui = ui, server = server)
